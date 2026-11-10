/*---------------------------------------------------------------------------*\
  onlineStatsOF6 functionObject for OpenFOAM.org v6
  Computes time-weighted online statistics and writes means/fluctuation terms.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvMeshFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace functionObjects
{

class onlineStatsOF6
:
    public fvMeshFunctionObject
{
    // Controls
    scalar timeStart_;
    scalar timeEnd_;

    // Helpers: replace '.' with '_' for accumulator object names
    static word safeName(const word& s)
    {
        word out(s);
        forAll(out, i)
        {
            if (out[i] == '.') out[i] = '_';
        }
        return out;
    }

    static word accName(const word& base, const word& fld)
    {
        return base + "_" + safeName(fld);
    }

    // Get or create and store a volScalarField accumulator
    volScalarField& getOrCreateVolScalar
    (
        const word& objName,
        const dimensionSet& dims
    )
    {
        if (mesh_.foundObject<volScalarField>(objName))
        {
            return const_cast<volScalarField&>
            (
                mesh_.lookupObject<volScalarField>(objName)
            );
        }

        volScalarField* fPtr = new volScalarField
        (
            IOobject
            (
                objName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dims, 0.0)
        );

        mesh_.objectRegistry::store(fPtr);
        return *fPtr;
    }

    volVectorField& getOrCreateVolVector
    (
        const word& objName,
        const dimensionSet& dims
    )
    {
        if (mesh_.foundObject<volVectorField>(objName))
        {
            return const_cast<volVectorField&>
            (
                mesh_.lookupObject<volVectorField>(objName)
            );
        }

        volVectorField* fPtr = new volVectorField
        (
            IOobject
            (
                objName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dims, vector::zero)
        );

        mesh_.objectRegistry::store(fPtr);
        return *fPtr;
    }

    volSymmTensorField& getOrCreateVolSymmTensor
    (
        const word& objName,
        const dimensionSet& dims
    )
    {
        if (mesh_.foundObject<volSymmTensorField>(objName))
        {
            return const_cast<volSymmTensorField&>
            (
                mesh_.lookupObject<volSymmTensorField>(objName)
            );
        }

        volSymmTensorField* fPtr = new volSymmTensorField
        (
            IOobject
            (
                objName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedSymmTensor("zero", dims, symmTensor::zero)
        );

        mesh_.objectRegistry::store(fPtr);
        return *fPtr;
    }

    volTensorField& getOrCreateVolTensor
    (
        const word& objName,
        const dimensionSet& dims
    )
    {
        if (mesh_.foundObject<volTensorField>(objName))
        {
            return const_cast<volTensorField&>
            (
                mesh_.lookupObject<volTensorField>(objName)
            );
        }

        volTensorField* fPtr = new volTensorField
        (
            IOobject
            (
                objName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedTensor("zero", dims, tensor::zero)
        );

        mesh_.objectRegistry::store(fPtr);
        return *fPtr;
    }

    // Accumulate scalar stats + scalar flux with U
    void accumulateScalarWithU
    (
        const word& fldName,
        const volVectorField& U,
        const volScalarField& phi,
        const scalar dt
    )
    {
        const word sumN  = accName("sum",  fldName);
        const word sum2N = accName("sum2", fldName);
        const word sumUN = accName("sumU", fldName);

        volScalarField& sumPhi  = getOrCreateVolScalar(sumN,  dimTime*phi.dimensions());
        volScalarField& sumPhi2 = getOrCreateVolScalar(sum2N, dimTime*sqr(phi.dimensions()));
        volVectorField& sumUPhi = getOrCreateVolVector(sumUN, dimTime*(U.dimensions()*phi.dimensions()));

        sumPhi.primitiveFieldRef()  += dt * phi.primitiveField();
        sumPhi2.primitiveFieldRef() += dt * sqr(phi.primitiveField());
        sumUPhi.primitiveFieldRef() += dt * (U.primitiveField() * phi.primitiveField());

        sumPhi.correctBoundaryConditions();
        sumPhi2.correctBoundaryConditions();
        sumUPhi.correctBoundaryConditions();
    }

    void writeScalarOutputs
    (
        const word& fldName,
        const volVectorField& UMean,
        volScalarField& sumDt
    )
    {
        const word sumN  = accName("sum",  fldName);
        const word sum2N = accName("sum2", fldName);
        const word sumUN = accName("sumU", fldName);

        if (!mesh_.foundObject<volScalarField>(sumN)) return;

        const volScalarField& sumPhi  = mesh_.lookupObject<volScalarField>(sumN);
        const volScalarField& sumPhi2 = mesh_.lookupObject<volScalarField>(sum2N);
        const volVectorField& sumUPhi = mesh_.lookupObject<volVectorField>(sumUN);

        volScalarField phiMean
        (
            IOobject(fldName + "Mean", mesh_.time().timeName(), mesh_,
                     IOobject::NO_READ, IOobject::AUTO_WRITE),
            sumPhi/sumDt
        );

        volScalarField phiPrime2Mean
        (
            IOobject(fldName + "Prime2Mean", mesh_.time().timeName(), mesh_,
                     IOobject::NO_READ, IOobject::AUTO_WRITE),
            (sumPhi2/sumDt) - sqr(phiMean)
        );

        volVectorField phiPrimeUPrimeMean
        (
            IOobject(fldName + "PrimeUPrimeMean", mesh_.time().timeName(), mesh_,
                     IOobject::NO_READ, IOobject::AUTO_WRITE),
            (sumUPhi/sumDt) - (UMean * phiMean)
        );

        phiMean.write();
        phiPrime2Mean.write();
        phiPrimeUPrimeMean.write();

        // Also write accumulators for restart
        const_cast<volScalarField&>(sumPhi).write();
        const_cast<volScalarField&>(sumPhi2).write();
        const_cast<volVectorField&>(sumUPhi).write();
    }

public:

    TypeName("onlineStatsOF6");

    onlineStatsOF6
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    )
    :
        fvMeshFunctionObject(name, runTime, dict),
        timeStart_(dict.lookupOrDefault<scalar>("timeStart", -GREAT)),
        timeEnd_(dict.lookupOrDefault<scalar>("timeEnd",   GREAT))
    {
        read(dict);
    }

    virtual bool read(const dictionary& dict)
    {
        fvMeshFunctionObject::read(dict);

        timeStart_ = dict.lookupOrDefault<scalar>("timeStart", -GREAT);
        timeEnd_   = dict.lookupOrDefault<scalar>("timeEnd",    GREAT);

        return true;
    }

    virtual bool execute()
    {
        const scalar t = mesh_.time().value();
        if (t < timeStart_ || t > timeEnd_) return true;

        const scalar dt = mesh_.time().deltaTValue();

        // Input fields
        if (!mesh_.foundObject<volVectorField>("U")) return true;
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

        // sumDt (stored as volScalarField for restart)
        volScalarField& sumDt = getOrCreateVolScalar("sumDt", dimTime);
        sumDt.primitiveFieldRef() += dt;
        sumDt.correctBoundaryConditions();

        // Velocity accumulators
        const word sumUName  = accName("sum", "U");
        const word sumUUName = accName("sum", "UU");

        volVectorField& sumU =
            getOrCreateVolVector(sumUName, dimTime*U.dimensions());

        volSymmTensorField& sumUU =
            getOrCreateVolSymmTensor(sumUUName, dimTime*sqr(U.dimensions()));

        sumU.primitiveFieldRef()  += dt * U.primitiveField();
        sumUU.primitiveFieldRef() += dt * symm(U.primitiveField()*U.primitiveField());

        sumU.correctBoundaryConditions();
        sumUU.correctBoundaryConditions();

        // Scalars
        if (mesh_.foundObject<volScalarField>("T"))
        {
            accumulateScalarWithU("T", U, mesh_.lookupObject<volScalarField>("T"), dt);
        }
        if (mesh_.foundObject<volScalarField>("p_rgh"))
        {
            accumulateScalarWithU("p_rgh", U, mesh_.lookupObject<volScalarField>("p_rgh"), dt);
        }
        if (mesh_.foundObject<volScalarField>("uMAGSqr"))
        {
            accumulateScalarWithU("uMAGSqr", U, mesh_.lookupObject<volScalarField>("uMAGSqr"), dt);
        }

        // Cross terms: U and force.disk_11
        if (mesh_.foundObject<volVectorField>("force.disk_11"))
        {
            const volVectorField& F = mesh_.lookupObject<volVectorField>("force.disk_11");

            const word sumFName  = accName("sum",   "force.disk_11");
            const word sumUFName = accName("sumUF", "force.disk_11");

            volVectorField& sumF =
                getOrCreateVolVector(sumFName, dimTime*F.dimensions());

            volTensorField& sumUF =
                getOrCreateVolTensor(sumUFName, dimTime*(U.dimensions()*F.dimensions()));

            sumF.primitiveFieldRef()  += dt * F.primitiveField();
            sumUF.primitiveFieldRef() += dt * (U.primitiveField()*F.primitiveField());

            sumF.correctBoundaryConditions();
            sumUF.correctBoundaryConditions();
        }

        return true;
    }

    virtual bool write()
    {
        // Need U and accumulators
        if (!mesh_.foundObject<volVectorField>("U")) return true;
        const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

        if (!mesh_.foundObject<volScalarField>("sumDt")) return true;
        volScalarField& sumDt = const_cast<volScalarField&>(mesh_.lookupObject<volScalarField>("sumDt"));

        const scalar denom = gMax(sumDt);
        if (denom <= SMALL) return true;

        const word sumUName  = accName("sum", "U");
        const word sumUUName = accName("sum", "UU");

        if (!mesh_.foundObject<volVectorField>(sumUName))  return true;
        if (!mesh_.foundObject<volSymmTensorField>(sumUUName)) return true;

        const volVectorField& sumU  = mesh_.lookupObject<volVectorField>(sumUName);
        const volSymmTensorField& sumUU = mesh_.lookupObject<volSymmTensorField>(sumUUName);

        // Means
        volVectorField UMean
        (
            IOobject("UMean", mesh_.time().timeName(), mesh_,
                     IOobject::NO_READ, IOobject::AUTO_WRITE),
            sumU/sumDt
        );
        UMean.write();

        volSymmTensorField UPrimeUPrimeMean
        (
            IOobject("UPrimeUPrimeMean", mesh_.time().timeName(), mesh_,
                     IOobject::NO_READ, IOobject::AUTO_WRITE),
            (sumUU/sumDt) - symm(UMean*UMean)
        );
        UPrimeUPrimeMean.write();

        // Scalars
        writeScalarOutputs("T",      UMean, sumDt);
        writeScalarOutputs("p_rgh",  UMean, sumDt);
        writeScalarOutputs("uMAGSqr",UMean, sumDt);

        // U-force cross terms: <U' outer F'>
        if (mesh_.foundObject<volVectorField>("force.disk_11"))
        {
            const word sumFName  = accName("sum",   "force.disk_11");
            const word sumUFName = accName("sumUF", "force.disk_11");

            if (mesh_.foundObject<volVectorField>(sumFName) &&
                mesh_.foundObject<volTensorField>(sumUFName))
            {
                const volVectorField& sumF  = mesh_.lookupObject<volVectorField>(sumFName);
                const volTensorField& sumUF = mesh_.lookupObject<volTensorField>(sumUFName);

                volVectorField FMean
                (
                    IOobject("force.disk_11Mean", mesh_.time().timeName(), mesh_,
                             IOobject::NO_READ, IOobject::AUTO_WRITE),
                    sumF/sumDt
                );
                FMean.write();

                volTensorField UForce_disk_11PrimeMean
                (
                    IOobject("UForce_disk_11PrimeMean", mesh_.time().timeName(), mesh_,
                             IOobject::NO_READ, IOobject::AUTO_WRITE),
                    (sumUF/sumDt) - (UMean*FMean)
                );
                UForce_disk_11PrimeMean.write();

                const_cast<volVectorField&>(sumF).write();
                const_cast<volTensorField&>(sumUF).write();
            }
        }

        // Write main accumulators for restart
        sumDt.write();
        const_cast<volVectorField&>(sumU).write();
        const_cast<volSymmTensorField&>(sumUU).write();

        return true;
    }
}; // end class onlineStatsOF6

defineTypeNameAndDebug(onlineStatsOF6, 0);
addToRunTimeSelectionTable(fvMeshFunctionObject, onlineStatsOF6, dictionary);

} // namespace functionObjects
} // namespace Foam