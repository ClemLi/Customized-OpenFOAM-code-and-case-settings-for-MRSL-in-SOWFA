/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Class
    Foam::fv::flyingActuationDiskSource

Description
    Actuation disk source for modeling multi-rotor system with lifting devices (MRSL)
    Made to be compatible with OpenFOAM v2506 and OpenFOAM-6

    Developed by YuanTso Li in Delft University of Technology, October, 2025
    Y.Li-18@tudelft.nl

    During the development, it is found that label and labelList is kind of "volatile"
    when running the code. Therefore, extra care should be paid when modifying the codes.
    Also, seems that if some of the variables will be contminated when being passed to funtions, 
    the reasons behind is unclear.

    On should avoid letting the actuator elements fall "exactly" at the boundary of cells. This
    may cause segmentation fault when excecuting a decomposed case.

SourceFiles
    flyingActuationDiskSource.C
    flyingActuationDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef flyingActuationDiskSource_H
#define flyingActuationDiskSource_H

#include "cellSetOption.H"
#include "Function1.H"
#include "volFieldsFwd.H"
#include "OFstream.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "labelList.H"
#include "scalarList.H"
#include "vectorList.H"
#include "interpolationCellPoint.H"
#include "Pstream.H"
#include "PstreamReduceOps.H" 


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class flyingActuationDiskSource Declaration
\*---------------------------------------------------------------------------*/

class flyingActuationDiskSource
:
    public cellSetOption
{

protected:

    // Protected data

    // Fields realted to IO and checker

        //- Output file stream
        OFstream* outputFileDetailBlade_;

        OFstream* outputFileDetailRotor_;

        OFstream* outputFileGeneral_;

        scalar timeStepToCheck_;

        scalar currentTimeStep_;

        label rotorForceActiveNumChecker_;

        label bladeForceActiveNumChecker_;

        label bladeExertingActiveListChecker_;

        label rotorExertingActiveListChecker_;

        scalar ringVelocityGlobalChecker_;

        scalar ringVelocityLocalChecker_;

        scalar currentProcsID_;

        scalar outPutTimeStart_;

        scalarListList outputFileDataGeneral_;

        scalarListList outputFileRotorInflowVelX_;
        scalarListList outputFileRotorInflowVelY_;
        scalarListList outputFileRotorInflowVelZ_;

        scalarListList outputFileRotorALForceX_;
        scalarListList outputFileRotorALForceY_;
        scalarListList outputFileRotorALForceZ_;

        scalarListList outputFileBladeInflowVelX_;
        scalarListList outputFileBladeInflowVelY_;
        scalarListList outputFileBladeInflowVelZ_;

        scalarListList outputFileBladeInflowVelOnPointX_;
        scalarListList outputFileBladeInflowVelOnPointY_;
        scalarListList outputFileBladeInflowVelOnPointZ_;

        scalarListList outputFileBladeALForceX_;
        scalarListList outputFileBladeALForceY_;
        scalarListList outputFileBladeALForceZ_;

        scalarListList outputFileBladeCl_;

        scalarListList outputFileBladeCd_;

        scalarListList outputFileBladePitch_;

        scalarList outputFileTime_;

        // controls the dumping interval of output files and the reset of forceFields
        label writeInterval_;

        label writeIntervalTrack_;


    // Fields related to MRSL chracteristics

        scalar rotorCtStar_;

        //- Actuator disk planar surface area [m2]
        scalar diskArea_;

        scalar widthOfDisk_;

        scalar heightOfDisk_;

        //- Surface-normal vector of the actuator disk for rotor
        vector diskDir_;

        vector freeStreamDir_;

        label bladeNumber_;

        scalar bladeChord_;

        //- Span direction of blades
        vector lateralDirection_;

        scalarList bladeTwist_;

    // Fields related to actuator elements chracteristics

        scalarList rotorFrontalAreaPerCell_;

        scalar bladeSpanPerAL_;

        vectorList rotorActuatorPointsPositions_;

        vectorList bladeActuatorPointsPositions_;


    // Fields related to forces. 
        // Lacal: local in each processor. 
        // Total: summing across all processors.
        // Without calc: considering the exacts AL points.
        // With calc: considering the smeared body force fields.

        //- Force field (per unit density per volume)
        volVectorField forceField_;

        vector totalForce_;

        vector localForce_;

        vector totalForceRotor_;

        vector localForceRotor_;

        vector totalForceBlade_;

        vector localForceBlade_;

        vector totalForceRotorCalc_;

        vector localForceRotorCalc_;

        vector totalForceBladeCalc_;

        vector localForceBladeCalc_;

        scalar rotorPowerCalc_;

        scalar bladePowerCalc_;


    // Fields related to sampling and smearing

        scalarList rotorSmearingLengthX_;
        scalarList rotorSmearingLengthY_;
        scalarList rotorSmearingLengthZ_;

        scalarList bladeSmearingLength_;

        labelList rotorCellID_;

        labelList bladeCellID_;

        vector inflowVelocity_;      //temporary storage for the blade element

        vector inflowVelocityOnPoint_;

        scalar velSamplingRadiusCoeff_;

        label velSamplingNumber_;

        scalarListList alphaVSClCd_;

        scalarList alphaList_;

        scalarList clList_;

        scalarList cdList_; 

        scalarList bladeAOATarget_;          // length = num of blade

        scalarList bladePitch_;              // length = num of AL points

        scalarList bladeAOA_;                // length = num of AL points

        scalarList bladeFlowAngle_;

        scalarList bladeAOAMidofThatBlade_;  // length = num of blade

        labelList bladeNumIndicator_;        // length = num of AL points

        labelList indexForBladeMid_;         // length = num of blade

        scalar timeStartAdjustPitch_;        // if value smaller than zero -> deactivated

        scalar bladePitchTol_;

        scalar bladePitchTimeScale_;

        scalar bladePitchSpeed_;             // unit: deg/sec

        scalarList misTargetAOAFiltered_;


    // Fields related parallelization

        // To accelerate the searching process with the pre-selecting feature
        labelListList rotorLabelListtoSmear_;

        labelListList bladeLabelListtoSmear_;

        // To accelerate the process by pre-calculating the smearing factor
        scalarListList rotorSmearFactorListList_;

        scalarListList bladeSmearFactorListList_;

        label bladeForceExertingActiveNum_;

        label rotorForceExertingActiveNum_;

        label bladeForceRecievingActiveNum_;

        label rotorForceRecievingActiveNum_;

        vectorList bladeALForceList_;

        vectorList rotorALForceList_;

        scalarList bladeALWhichProcs_;

        scalarList rotorALWhichProcs_;

        labelList bladeExertingActiveList_;

        labelList rotorExertingActiveList_;

        scalarList bladeClList_;                      // length = number of blade AL points

        scalarList bladeCdList_;                      // length = number of blade AL points

        vectorList inflowVelocityRotorList_;          // length = number of rotro AL points

        vectorList inflowVelocityList_;               // length = number of blade AL points

        vectorList inflowVelocityOnPointList_;        // length = number of blade AL points

        labelListList inflowVelocityRingSampleList_;  // length = number of blade AL points * numbeer of ring velocity sampling number



    // Member Functions

        Foam::label findCell
        (
            const point& location
        );

        //- Create the performance output file
        virtual void createOutputFile();

        //- Interpolate a scalar value
        scalar interpolate
        (
            scalar xNew,
            scalarList& xOld,
            scalarList& yOld
        );

private:

    // Private Member Functions

        //- Compute variable-scaling actuator disk method
        template<class AlphaFieldType, class RhoFieldType>
        void makeWindTurbineWakesFly
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn,
            interpolationCellPoint<vector> UInterp
        );

        //- Disallow default bitwise copy construct
        flyingActuationDiskSource(const flyingActuationDiskSource&);

        //- Disallow default bitwise assignment
        void operator=(const flyingActuationDiskSource&);


public:


    //- Pre-selecting cells to smear from the start instead of calculating them every time
            void selectingCelltoSmear(
                vector centerToSmear,
                label currentIndex,
                scalar smearingLengthX,
                scalar smearingLengthY,
                scalar smearingLengthZ,
                label rotorOrBlade      // 0 stands for rotor, 1 stands for blade
            );


    //- Apply force field based on force vector
            void applyForceField(fvMatrix<vector>& eqn, const vectorField& U, vector cellCenterNow, vector forceVectorToApply, 
                labelList domainCells, scalarList smearFactorForDomainCells, label rotorOrBlade);


    //- Runtime type information
    TypeName("flyingActuationDiskSource");


    // Constructors

        //- Construct from components
        flyingActuationDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );



    //- Destructor
    virtual ~flyingActuationDiskSource();


    // Member Functions

        // Access

            //- Normal disk direction
            const vector& diskDir() const
            {
                return diskDir_;
            }

            //- Disk area
            scalar diskArea() const
            {
                return diskArea_;
            }


         // Add explicit and implicit contributions

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "flyingActuationDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
