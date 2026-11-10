/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

\*---------------------------------------------------------------------------*/

#include "flyingActuationDiskSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(flyingActuationDiskSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        flyingActuationDiskSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::flyingActuationDiskSource::~flyingActuationDiskSource()
{}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::flyingActuationDiskSource::flyingActuationDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    diskDir_(coeffs_.lookup("diskDir")),
    diskArea_(readScalar(coeffs_.lookup("diskArea"))),
    rotorCtStar_(readScalar(coeffs_.lookup("rotorCtStar"))),
    widthOfDisk_(readScalar(coeffs_.lookup("widthOfDisk"))),
    heightOfDisk_(readScalar(coeffs_.lookup("heightOfDisk"))),
    totalForce_(Zero),
    localForce_(Zero),
    totalForceRotor_(Zero),
    localForceRotor_(Zero),
    totalForceBlade_(Zero),
    localForceBlade_(Zero),
    bladeChord_(readScalar(coeffs_.lookup("bladeChord"))),
    bladeNumber_(readLabel(coeffs_.lookup("bladeNumber"))),
    bladeSpanPerAL_(readScalar(coeffs_.lookup("bladeSpanPerAL"))),
    velSamplingRadiusCoeff_(readScalar(coeffs_.lookup("velSamplingRadiusCoeff"))),
    velSamplingNumber_(readLabel(coeffs_.lookup("velSamplingNumber"))),
    timeStartAdjustPitch_(readScalar(coeffs_.lookup("timeStartAdjustPitch"))),
    bladePitchTol_(readScalar(coeffs_.lookup("bladePitchTol"))),
    inflowVelocity_(Zero),
    bladePitchSpeed_(readScalar(coeffs_.lookup("bladePitchSpeed"))),
    bladePitchTimeScale_(readScalar(coeffs_.lookup("bladePitchTimeScale"))),
    timeStepToCheck_(readScalar(coeffs_.lookup("timeStepToCheck"))),
    outPutTimeStart_(readScalar(coeffs_.lookup("outPutTimeStart"))),
    writeInterval_(readLabel(coeffs_.lookup("writeInterval"))),
    forceField_
    (
        IOobject
        (
            "force." + name_,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "force",
            dimForce/dimVolume/dimDensity,
            vector::zero
        )
    )
{

    currentProcsID_ = round( Pstream::myProcNo() );

    coeffs_.lookup("bladeActuatorPointsPositions") >> bladeActuatorPointsPositions_;

    coeffs_.lookup("rotorActuatorPointsPositions") >> rotorActuatorPointsPositions_;

    coeffs_.lookup("alphaVSClCd") >> alphaVSClCd_;

    alphaList_ = scalarList(alphaVSClCd_.size(), 0.0);

    clList_ = scalarList(alphaVSClCd_.size(), 0.0);

    cdList_ = scalarList(alphaVSClCd_.size(), 0.0);

    misTargetAOAFiltered_ = scalarList(bladeNumber_, 0.0);

    forAll(alphaVSClCd_, jj)
    {
        alphaList_[jj] = alphaVSClCd_[jj][0];
        clList_[jj] = alphaVSClCd_[jj][1];
        cdList_[jj] = alphaVSClCd_[jj][2];
    }

    dict.readIfPresent("bladeSmearingLength", bladeSmearingLength_);

    dict.readIfPresent("rotorSmearingLengthX", rotorSmearingLengthX_);
    dict.readIfPresent("rotorSmearingLengthY", rotorSmearingLengthY_);
    dict.readIfPresent("rotorSmearingLengthZ", rotorSmearingLengthZ_);

    dict.readIfPresent("rotorFrontalAreaPerCell", rotorFrontalAreaPerCell_);

    dict.readIfPresent("bladeTwist", bladeTwist_);

    dict.readIfPresent("bladePitch", bladePitch_);

    dict.readIfPresent("bladeAOATarget", bladeAOATarget_);

    dict.readIfPresent("bladeNumIndicator", bladeNumIndicator_);

    dict.readIfPresent("indexForBladeMid", indexForBladeMid_);

    currentTimeStep_ = 0.0;

    coeffs_.lookup("fields") >> fieldNames_;
    
    applied_.setSize(fieldNames_.size(), false);

    Info<< "    - creating actuation MRSL zone: " << this->name() << endl;

    forceField_.write();

    vector freeStreamDirTemp;
    dict.readIfPresent("freeStreamDir", freeStreamDirTemp);
    freeStreamDir_ = freeStreamDirTemp / mag(freeStreamDirTemp);
    if (mag(freeStreamDir_) < VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Free stream direction vector is zero: "
            << "freeStreamDir = " << freeStreamDir_
            << exit(FatalIOError);
    }

    vector lateralDirectionTemp;
    dict.readIfPresent("lateralDirection", lateralDirectionTemp);
    lateralDirection_ = lateralDirectionTemp / mag(lateralDirectionTemp);
    if (mag(lateralDirection_) < VSMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Free stream direction vector is zero: "
            << "lateralDirection = " << lateralDirection_
            << exit(FatalIOError);
    }


    /////////////////////// Parallelization preparation ///////////////////////


    rotorLabelListtoSmear_.setSize( rotorActuatorPointsPositions_.size() );
    forAll(rotorLabelListtoSmear_, idxFilling)
    {
        rotorLabelListtoSmear_[idxFilling] = Foam::labelList(1, Foam::label(-1) );
    }

    rotorSmearFactorListList_.setSize( rotorActuatorPointsPositions_.size() );
    forAll(rotorSmearFactorListList_, idxFilling)
    {
        rotorSmearFactorListList_[idxFilling] = Foam::scalarList(1, Foam::scalar(0.0) );
    }

    bladeLabelListtoSmear_.setSize( bladeActuatorPointsPositions_.size() );
    forAll(bladeLabelListtoSmear_, idxFilling)
    {
        bladeLabelListtoSmear_[idxFilling] = Foam::labelList(1, Foam::label(-1) );
    }

    bladeSmearFactorListList_.setSize( bladeActuatorPointsPositions_.size() );
    forAll(bladeSmearFactorListList_, idxFilling)
    {
        bladeSmearFactorListList_[idxFilling] = Foam::scalarList(1, Foam::scalar(0.0) );
    }


    inflowVelocityRingSampleList_.setSize( bladeActuatorPointsPositions_.size() );
    forAll(inflowVelocityRingSampleList_, idxFilling)
    {
        inflowVelocityRingSampleList_[idxFilling] = Foam::labelList(velSamplingNumber_, Foam::label(-1) );
    }

    rotorExertingActiveList_ = labelList(rotorActuatorPointsPositions_.size(), 0);
    bladeExertingActiveList_ = labelList(bladeActuatorPointsPositions_.size(), 0);

    rotorForceExertingActiveNum_ = 0;
    bladeForceExertingActiveNum_ = 0;
    rotorForceRecievingActiveNum_ = 0;
    bladeForceRecievingActiveNum_ = 0;

    
    rotorALWhichProcs_ = scalarList(rotorActuatorPointsPositions_.size(), 0.0);
    bladeALWhichProcs_ = scalarList(bladeActuatorPointsPositions_.size(), 0.0);

    bladeAOA_ = scalarList(bladeActuatorPointsPositions_.size(), 0.0);
    bladeAOAMidofThatBlade_ = scalarList(bladeNumber_, 0.0);

    rotorCellID_ = labelList(rotorActuatorPointsPositions_.size(), -1);
    bladeCellID_ = labelList(bladeActuatorPointsPositions_.size(), -1);

    /////////////////////////////////// Parallelization for rotor ///////////////////////////////////

    /*
    Identifying which processors to be applied the force and pre-selecting cells to apply the body force
    fields of the rotor (thrust)
    */

    label rotorCurrentIndex = 0;

    for ( 
        label rotorActuatorPointsCounter = 0; 
        rotorActuatorPointsCounter < rotorActuatorPointsPositions_.size(); 
        rotorActuatorPointsCounter++
    )
    {

        vector centerToSmear = rotorActuatorPointsPositions_[rotorActuatorPointsCounter];
        label rotorOrBlade = 0;

        selectingCelltoSmear(centerToSmear, rotorCurrentIndex, rotorSmearingLengthX_[rotorActuatorPointsCounter],
            rotorSmearingLengthY_[rotorActuatorPointsCounter], rotorSmearingLengthZ_[rotorActuatorPointsCounter], rotorOrBlade);

        rotorCurrentIndex += 1;
    }

    /*
    Identifying which processors that exert the force of the rotor (thrust)
    */

    
    for (
        label rotorActuatorPointsCounter = 0; 
        rotorActuatorPointsCounter < rotorActuatorPointsPositions_.size(); 
        rotorActuatorPointsCounter++
    )
    {
        vector vectorHolder = rotorActuatorPointsPositions_[rotorActuatorPointsCounter];
        label labelHolder = findCell(vectorHolder);

        if ( labelHolder > -1 )
        {
            rotorForceExertingActiveNum_ += 1;
            rotorExertingActiveList_[rotorActuatorPointsCounter] = 1;
            rotorALWhichProcs_[rotorActuatorPointsCounter] = currentProcsID_;

            // Pre-selecting
            rotorCellID_[rotorActuatorPointsCounter] = labelHolder;

            Pout << name_ << ": rotorForceExerting element " 
                << rotorActuatorPointsCounter << " at " << vectorHolder << endl;
        }

    }

    for 
    (
        label actuatorPointIndex = 0; 
        actuatorPointIndex < rotorActuatorPointsPositions_.size();
        actuatorPointIndex++
    )
    {
        reduce( rotorALWhichProcs_[actuatorPointIndex], sumOp<scalar>() );
    }

    Pout << name_ << ": rotorForceExertingActiveNum: " << rotorForceExertingActiveNum_ << endl;
    Pout << name_ << ": rotorForceRecievingActiveNum: " << rotorForceRecievingActiveNum_ << endl;

    /////////////////////////////////// Parallelization for blades ///////////////////////////////////

    /*
    Identifying which processors to be applied the force and pre-selecting cells to apply the body force
    fields of the blades (lift and drag)
    */

    label bladeCurrentIndex = 0;

    for ( 
        label bladeActuatorPointsCounter = 0; 
        bladeActuatorPointsCounter < bladeActuatorPointsPositions_.size(); 
        bladeActuatorPointsCounter++
    )
    {
        vector centerToSmear = bladeActuatorPointsPositions_[bladeActuatorPointsCounter];
        label rotorOrBlade = 1;

        selectingCelltoSmear(centerToSmear, bladeCurrentIndex, bladeSmearingLength_[bladeActuatorPointsCounter],
            bladeSmearingLength_[bladeActuatorPointsCounter], bladeSmearingLength_[bladeActuatorPointsCounter], rotorOrBlade);

        bladeCurrentIndex += 1;

    }

    /*
    Identifying which processors that exert the force of the blades (lift and drags)
    */

    for (
        label bladeActuatorPointsCounter = 0; 
        bladeActuatorPointsCounter < bladeActuatorPointsPositions_.size(); 
        bladeActuatorPointsCounter++
    )
    {
        vector vectorHolder = bladeActuatorPointsPositions_[bladeActuatorPointsCounter];
        label labelHolder = findCell(vectorHolder);

        if ( labelHolder > -1 )
        {
            bladeForceExertingActiveNum_ += 1;
            bladeExertingActiveList_[bladeActuatorPointsCounter] = 1;
            bladeALWhichProcs_[bladeActuatorPointsCounter] = currentProcsID_;

            // Pre-selecting
            bladeCellID_[bladeActuatorPointsCounter] = labelHolder;

            Pout << name_ << ": bladeForceExerting element " 
                << bladeActuatorPointsCounter << " at " << vectorHolder << endl;
        }
    }


    for 
    (
        label actuatorPointIndex = 0; 
        actuatorPointIndex < bladeActuatorPointsPositions_.size();
        actuatorPointIndex++
    )
    {
        reduce( bladeALWhichProcs_[actuatorPointIndex], sumOp<scalar>() );
    }

    /*
    Identifying which processors to select the velocity fields for the blades if using ring-like line sampling.
    Also, the cells are pre-selected.
    */

    ringVelocityGlobalChecker_ = 0.0;

    for
    (
        label actuatorPointIndex = 0; 
        actuatorPointIndex < bladeActuatorPointsPositions_.size();
        actuatorPointIndex++
    )
    {
        ringVelocityLocalChecker_ = 0.0;

        for( label velSampleIndex = 0; velSampleIndex < velSamplingNumber_; velSampleIndex++ )
        {

            scalar velSamplingRadius = bladeSmearingLength_[actuatorPointIndex] * velSamplingRadiusCoeff_;

            // distribute the points evenly in terms of angular distance
            scalar pointAngle = Foam::constant::mathematical::pi * 2.0 * velSampleIndex / velSamplingNumber_;
            scalar freeStreamDist = velSamplingRadius * Foam::cos(pointAngle);
            scalar verticalDist = velSamplingRadius * Foam::sin(pointAngle);
            vector samplePoint = bladeActuatorPointsPositions_[actuatorPointIndex]
                + freeStreamDist * freeStreamDir_ + verticalDist * (freeStreamDir_ ^ lateralDirection_);

            inflowVelocityRingSampleList_[actuatorPointIndex][velSampleIndex] = findCell(samplePoint);

            if ( inflowVelocityRingSampleList_[actuatorPointIndex][velSampleIndex] > -1 )
            {
                ringVelocityLocalChecker_ = ringVelocityLocalChecker_ + 1.0;
                ringVelocityGlobalChecker_ = ringVelocityGlobalChecker_ + 1.0;
            }
        }

        // checker
        scalar ringVelocityLocalCheckerDeclared = returnReduce( ringVelocityLocalChecker_, sumOp<scalar>() );

        if ( round( ringVelocityLocalCheckerDeclared ) != velSamplingNumber_ )
        {
            FatalErrorIn("velocity sampling points for blade AL point fails: local")
                << name_ << ": reduced ringVelocityLocalChecker_ at point " << actuatorPointIndex
                << " did not pass, returnReduce( ringVelocityLocalChecker_, sumOp<scalar>() ) is " 
                << ringVelocityLocalCheckerDeclared << ", but it should be "
                << velSamplingNumber_ << endl
                << "Note that the actuator elements (velocity sampling points as well) should not be right at the border of a cell, "
                << "this may cause issues as the program will have difficuties in decinging which processor is holding the actuator elements"
                << endl
                << abort(FatalError);
        }
        else if ( Pstream::master() )
        {
            Pout << name_ << ": reduced ringVelocityLocalChecker_ at point " << actuatorPointIndex
                << " passed, returnReduce( ringVelocityLocalChecker_, sumOp<scalar>() ) is " 
                << ringVelocityLocalCheckerDeclared << ", as it should be." << endl;
        }

    }

    Pout << name_ << ": bladeForceExertingActiveNum: " << bladeForceExertingActiveNum_ << endl;
    Pout << name_ << ": bladeForceRecievingActiveNum: " << bladeForceRecievingActiveNum_ << endl;

    /////////////////////////////////// Create output files ///////////////////////////////////

    createOutputFile();

    /////////////////////////////////// Some fields initialization ///////////////////////////////////

    // for computations

    localForceRotor_ = vector(0.0, 0.0, 0.0);
    localForceRotorCalc_ = vector(0.0, 0.0, 0.0);
    localForceBlade_ = vector(0.0, 0.0, 0.0);
    localForceBladeCalc_ = vector(0.0, 0.0, 0.0);

    rotorALForceList_ = vectorList(rotorActuatorPointsPositions_.size(), vector::zero);
    bladeALForceList_ = vectorList(bladeActuatorPointsPositions_.size(), vector::zero);
    inflowVelocityRotorList_ = vectorList(rotorActuatorPointsPositions_.size(), vector::zero);

    bladeFlowAngle_ = scalarList(bladeActuatorPointsPositions_.size(), 0.0);
    bladeAOA_ = scalarList(bladeActuatorPointsPositions_.size(), 0.0);

    bladeClList_ = scalarList(bladeActuatorPointsPositions_.size(), 0.0);
    bladeCdList_ = scalarList(bladeActuatorPointsPositions_.size(), 0.0);

    inflowVelocityList_ = vectorList(bladeActuatorPointsPositions_.size(), vector::zero);
    inflowVelocityOnPointList_ = vectorList(bladeActuatorPointsPositions_.size(), vector::zero);

    // for outputting

    writeIntervalTrack_ = 0;

    // general

    outputFileTime_ = Foam::scalarList( writeInterval_, Foam::scalar(0.0) );

    outputFileDataGeneral_.setSize( writeInterval_ );

    for 
    (
        label idxFilling = 0; 
        idxFilling < writeInterval_;
        idxFilling++
    )
    {
        outputFileDataGeneral_[idxFilling] = Foam::scalarList( 17, Foam::scalar(0.0) );
    }

    // rotor

    outputFileRotorInflowVelX_.setSize( writeInterval_ );
    outputFileRotorInflowVelY_.setSize( writeInterval_ );
    outputFileRotorInflowVelZ_.setSize( writeInterval_ );

    outputFileRotorALForceX_.setSize( writeInterval_ );
    outputFileRotorALForceY_.setSize( writeInterval_ );
    outputFileRotorALForceZ_.setSize( writeInterval_ );

    for 
    (
        label idxFilling = 0; 
        idxFilling < writeInterval_;
        idxFilling++
    )
    {
        outputFileRotorInflowVelX_[idxFilling] = Foam::scalarList( rotorActuatorPointsPositions_.size(), 0.0 );
        outputFileRotorInflowVelY_[idxFilling] = Foam::scalarList( rotorActuatorPointsPositions_.size(), 0.0 );
        outputFileRotorInflowVelZ_[idxFilling] = Foam::scalarList( rotorActuatorPointsPositions_.size(), 0.0 );

        outputFileRotorALForceX_[idxFilling] = Foam::scalarList( rotorActuatorPointsPositions_.size(), 0.0 );
        outputFileRotorALForceY_[idxFilling] = Foam::scalarList( rotorActuatorPointsPositions_.size(), 0.0 );
        outputFileRotorALForceZ_[idxFilling] = Foam::scalarList( rotorActuatorPointsPositions_.size(), 0.0 );
    }

    // blade

    outputFileBladeInflowVelX_.setSize( writeInterval_ );
    outputFileBladeInflowVelY_.setSize( writeInterval_ );
    outputFileBladeInflowVelZ_.setSize( writeInterval_ );

    outputFileBladeInflowVelOnPointX_.setSize( writeInterval_ );
    outputFileBladeInflowVelOnPointY_.setSize( writeInterval_ );
    outputFileBladeInflowVelOnPointZ_.setSize( writeInterval_ );

    outputFileBladeALForceX_.setSize( writeInterval_ );
    outputFileBladeALForceY_.setSize( writeInterval_ );
    outputFileBladeALForceZ_.setSize( writeInterval_ );

    outputFileBladeCl_.setSize( writeInterval_ );
    outputFileBladeCd_.setSize( writeInterval_ );
    outputFileBladePitch_.setSize( writeInterval_ );

    for 
    (
        label idxFilling = 0; 
        idxFilling < writeInterval_;
        idxFilling++
    )
    {
        outputFileBladeInflowVelX_[idxFilling] = Foam::scalarList( bladeActuatorPointsPositions_.size(), 0.0 );
        outputFileBladeInflowVelY_[idxFilling] = Foam::scalarList( bladeActuatorPointsPositions_.size(), 0.0 );
        outputFileBladeInflowVelZ_[idxFilling] = Foam::scalarList( bladeActuatorPointsPositions_.size(), 0.0 );

        outputFileBladeInflowVelOnPointX_[idxFilling] = Foam::scalarList( bladeActuatorPointsPositions_.size(), 0.0 );
        outputFileBladeInflowVelOnPointY_[idxFilling] = Foam::scalarList( bladeActuatorPointsPositions_.size(), 0.0 );
        outputFileBladeInflowVelOnPointZ_[idxFilling] = Foam::scalarList( bladeActuatorPointsPositions_.size(), 0.0 );

        outputFileBladeALForceX_[idxFilling] = Foam::scalarList( bladeActuatorPointsPositions_.size(), 0.0 );
        outputFileBladeALForceY_[idxFilling] = Foam::scalarList( bladeActuatorPointsPositions_.size(), 0.0 );
        outputFileBladeALForceZ_[idxFilling] = Foam::scalarList( bladeActuatorPointsPositions_.size(), 0.0 );

        outputFileBladeCl_[idxFilling] = Foam::scalarList( bladeActuatorPointsPositions_.size(), Foam::scalar(0.0) );
        outputFileBladeCd_[idxFilling] = Foam::scalarList( bladeActuatorPointsPositions_.size(), Foam::scalar(0.0) );
        outputFileBladePitch_[idxFilling] = Foam::scalarList( bladeActuatorPointsPositions_.size(), Foam::scalar(0.0) );
    }


    /////////////////////////////////// Some checking ///////////////////////////////////

    scalar ringVelocityGlobalCheckerDeclared = returnReduce( ringVelocityGlobalChecker_, sumOp<scalar>() );

    if ( round( ringVelocityGlobalCheckerDeclared ) != velSamplingNumber_ * bladeActuatorPointsPositions_.size() )
    {
        FatalErrorIn("velocity sampling points for blade AL point fails: global")
            << name_ << ": reduced ringVelocityGlobalChecker_ did not pass, " <<
            "returnReduce( ringVelocityGlobalChecker_, sumOp<scalar>() ) is " 
            << ringVelocityGlobalCheckerDeclared << ", but it should be "
            << velSamplingNumber_ * bladeActuatorPointsPositions_.size() << endl
            << "Note that the actuator elements (velocity sampling points as well) should not be right at the border of a cell, "
            << "this may cause issues as the program will have difficuties in decinging which processor is holding the actuator elements"
            << endl
            << abort(FatalError);
    }
    else if ( Pstream::master() )
    {
        Pout << name_ << ": reduced ringVelocityGlobalChecker_ at point passed, "
            << "returnReduce( ringVelocityGlobalChecker_, sumOp<scalar>() ) is " 
            << ringVelocityGlobalCheckerDeclared << ", as it should be." << endl;
    }

    for 
    (
        label actuatorPointIndex = 0; 
        actuatorPointIndex < rotorActuatorPointsPositions_.size();
        actuatorPointIndex++
    )
    {
        if ( returnReduce( rotorExertingActiveList_[actuatorPointIndex], sumOp<label>() ) != 1)
        {
            FatalErrorIn("force exerting AL number of rotor detail")
                << "reduced rotorExertingActiveList_ of " << name_
                << " did not pass, returnReduce( rotorExertingActiveList_[" << actuatorPointIndex 
                << "], sumOp<label>() ) is " 
                << returnReduce( rotorExertingActiveList_[actuatorPointIndex], sumOp<label>() ) << endl
                << "Note that the actuator elements (velocity sampling points as well) should not be right at the border of a cell, "
                << "this may cause issues as the program will have difficuties in decinging which processor is holding the actuator elements"
                << endl
                << abort(FatalError);
        }
    }

    for 
    (
        label actuatorPointIndex = 0; 
        actuatorPointIndex < bladeActuatorPointsPositions_.size();
        actuatorPointIndex++
    )
    {
        if ( returnReduce( bladeExertingActiveList_[actuatorPointIndex], sumOp<label>() ) != 1)
        {
            FatalErrorIn("force exerting AL number of blade detail")
                << "reduced bladeExertingActiveList_ of " << name_
                << " did not pass, returnReduce( bladeExertingActiveList_[" << actuatorPointIndex 
                << "], sumOp<label>() ) is " 
                << returnReduce( bladeExertingActiveList_[actuatorPointIndex], sumOp<label>() ) << endl
                << "Note that the actuator elements (velocity sampling points as well) should not be right at the border of a cell, "
                << "this may cause issues as the program will have difficuties in decinging which processor is holding the actuator elements"
                << endl
                << abort(FatalError);
        }
    }

    rotorExertingActiveListChecker_ = 0;
    bladeExertingActiveListChecker_ = 0;


    for 
    (
        label actuatorPointIndex = 0; 
        actuatorPointIndex < rotorActuatorPointsPositions_.size();
        actuatorPointIndex++
    )
    {
        rotorExertingActiveListChecker_ = rotorExertingActiveListChecker_ + rotorExertingActiveList_[actuatorPointIndex];
    }


    for 
    (
        label actuatorPointIndex = 0; 
        actuatorPointIndex < bladeActuatorPointsPositions_.size();
        actuatorPointIndex++
    )
    {
        bladeExertingActiveListChecker_ = bladeExertingActiveListChecker_ + bladeExertingActiveList_[actuatorPointIndex];
    }

    reduce( rotorExertingActiveListChecker_, sumOp<label>() );
    reduce( bladeExertingActiveListChecker_, sumOp<label>() );


    if (Pstream::master())
    {

        if ( rotorExertingActiveListChecker_ !=  rotorActuatorPointsPositions_.size() )
        {
            FatalErrorIn("force exerting AL number of rotor")
                << "The total force exerting AL number of rotor found in program for " << name_
                << " is " << rotorExertingActiveListChecker_ << ", but it should be "
                << rotorActuatorPointsPositions_.size() << endl
                << "recommend to check if all the AL points are located in cell center, "
                << "otherwise the program don't know where the point lies in"
                << abort(FatalError);
        }
        else
        {
            Pout << "The total force exerting AL number of rotor found in program for " << name_
                << " is " << rotorExertingActiveListChecker_ << endl;
        }

        if ( bladeExertingActiveListChecker_ !=  bladeActuatorPointsPositions_.size() )
        {
            FatalErrorIn("force exerting AL number of blade")
                << "The total force exerting AL number of blade found in program for " << name_
                << " is " << bladeExertingActiveListChecker_ << ", but it should be "
                << bladeActuatorPointsPositions_.size() << endl
                << "recommend to check if all the AL points are located in cell center, "
                << "otherwise the program don't know where the point lies in"
                << abort(FatalError);
        }
        else
        {
            Pout << "The total force exerting AL number of blade found in program for " << name_
                << " is " << bladeExertingActiveListChecker_ << endl;
        }

    }

    Pout << endl;
    Pout << "howManyALPointConsidered (" << name_ << " Rotor): " << rotorCurrentIndex << endl;
    Pout << "howManyALPointConsidered (" << name_ << " Blade): " << bladeCurrentIndex << endl;

    /////////////////////////////////// Summary of the object ///////////////////////////////////


    Info << endl;
    Info << "Summary of " << name_ << ":" << endl;
    Info << endl;
    Info << "rotorCtStar: " << rotorCtStar_ << endl;
    Info << "widthOfDisk: " << widthOfDisk_ << endl;
    Info << "heightOfDisk: " << heightOfDisk_ << endl;
    Info << "freeStreamDir: " << freeStreamDir_ << endl;
    Info << "lateralDirection: " << lateralDirection_ << endl;
    Info << "rotorSmearingLengthX[0]: " << rotorSmearingLengthX_[0] << endl;
    Info << "rotorSmearingLengthY[0]: " << rotorSmearingLengthY_[0] << endl;
    Info << "rotorSmearingLengthZ[0]: " << rotorSmearingLengthZ_[0] << endl;
    Info << "diskDir: " << diskDir_ << endl;
    Info << "rotorFrontalAreaPerCell[0]: " << rotorFrontalAreaPerCell_[0] << endl;
    Info << endl;
    Info << endl;
    Info << "bladeSmearingLength[0]: " << bladeSmearingLength_[0] << endl;
    Info << "bladeNumber: " << bladeNumber_ << endl;
    Info << "bladeChord: " << bladeChord_ << endl;
    Info << "bladeSpanPerAL: " << bladeSpanPerAL_ << endl;
    Info << "bladeTwist_[0]: " << bladeTwist_[0] << " deg" << endl;
    Info << "velSamplingRadiusCoeff: " << velSamplingRadiusCoeff_ << endl;
    Info << "velSamplingNumber: " << velSamplingNumber_ << endl;

    Info << "clList_[4]: " << clList_[4] << endl;
    Info << "bladeAOATarget_[0]: " << bladeAOATarget_[0] << endl;
    Info << "timeStartAdjustPitch: " << timeStartAdjustPitch_ << endl;
    Info << "bladePitchTol: " << bladePitchTol_ << endl;
    Info << "bladePitchSpeed: " << bladePitchSpeed_ << endl;
    Info << "bladePitchTimeScale: " << bladePitchTimeScale_ << endl;

    Info << "timeStepToCheck_: " << timeStepToCheck_ << endl;
    Info << "writeInterval_: " << writeInterval_ << endl;

    Info << endl;
    Info << endl;

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::flyingActuationDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    currentTimeStep_ = currentTimeStep_ + 1.0;

    //Somehow tricky but interpolationCellPoint<vector> does not work inside if statement

    const volVectorField& Uin(eqn.psi());
    interpolationCellPoint<vector> UInterp(Uin);

    makeWindTurbineWakesFly(geometricOneField(), geometricOneField(), eqn, UInterp);
}


void Foam::fv::flyingActuationDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    NotImplemented;
}


bool Foam::fv::flyingActuationDiskSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        dict.readIfPresent("diskArea", diskArea_);
        if (diskArea_ < VSMALL)
        {
            FatalIOErrorInFunction(dict)
                << "Actuator disk has zero area: "
                << "diskArea = " << diskArea_
                << exit(FatalIOError);
        }

        vector diskDirTemp;
        dict.readIfPresent("diskDir", diskDirTemp);
        diskDir_ = diskDirTemp / mag(diskDirTemp);
        if (mag(diskDir_) < VSMALL)
        {
            FatalIOErrorInFunction(dict)
                << "Actuator disk surface-normal vector is zero: "
                << "diskDir = " << diskDir_
                << exit(FatalIOError);
        }

        Info << "Done read()." << endl;

        return true;
    }

    return false;
}


Foam::label Foam::fv::flyingActuationDiskSource::findCell
(
    const point& location
)
{
    return mesh_.findCell(location);
}


// This member function poses risk to mess up with rotorALWhichProcs_ for unknown reasons
void Foam::fv::flyingActuationDiskSource::selectingCelltoSmear(
    vector centerToSmear, label currentIndex, scalar smearingLengthX,
    scalar smearingLengthY, scalar smearingLengthZ, label rotorOrBlade)
{

    scalar largestSmearingLength = max(smearingLengthX, smearingLengthY);
    largestSmearingLength = max(largestSmearingLength, smearingLengthZ);

    scalar sphereRadius = largestSmearingLength * 3;
    
    scalar smearCellCounter = 0;

    labelList selectedCelltoSmear;

    scalarList smearFactorForTheSelectedCell;

    forAll(mesh_.cells(), cellI)
    {
        vector smearDisplacementVec = mesh_.C()[cellI] - centerToSmear;
        
        scalar dis = mag( smearDisplacementVec );

        scalar disX = smearDisplacementVec.x();
        scalar disY = smearDisplacementVec.y();
        scalar disZ = smearDisplacementVec.z();

        if(dis <= sphereRadius)
        {
            selectedCelltoSmear.append(cellI);
            smearCellCounter += 1;

            scalar factor = Foam::exp( -Foam::sqr(disX/smearingLengthX) -Foam::sqr(disY/smearingLengthY) -Foam::sqr(disZ/smearingLengthZ) )
                      / ( (smearingLengthX*smearingLengthY*smearingLengthZ)  // it is divided by Volume at here -> force becomes force density
                      * Foam::pow(Foam::constant::mathematical::pi, 1.5));

            smearFactorForTheSelectedCell.append(factor);
        }
    }


    if (smearCellCounter > 0)
    {
        if(rotorOrBlade == 0)
        {
            rotorLabelListtoSmear_[currentIndex] = selectedCelltoSmear;
            rotorSmearFactorListList_[currentIndex] = smearFactorForTheSelectedCell;
            rotorForceRecievingActiveNum_ = rotorForceRecievingActiveNum_ + 1;

            Pout << "select cell to smear, currentIndex (" << name_ << " rotor): " 
                << currentIndex << ", Position at: " << centerToSmear << endl;
            Pout << "smearCellCounter during initializing (" << name_ << " rotor): " << smearCellCounter << endl;


        }

        if(rotorOrBlade == 1)
        {
            bladeLabelListtoSmear_[currentIndex] = selectedCelltoSmear;
            bladeSmearFactorListList_[currentIndex] = smearFactorForTheSelectedCell;
            bladeForceRecievingActiveNum_ = bladeForceRecievingActiveNum_ + 1;

            Pout << "select cell to smear, currentIndex (" << name_ << " blade): " 
                << currentIndex << ", Position at: " << centerToSmear << endl;
            Pout << "smearCellCounter during initializing (" << name_ << " blade): " << smearCellCounter << endl;
        }
    }

}


void Foam::fv::flyingActuationDiskSource::createOutputFile()
{

    // blade detail

    fileName dirDetailBlade;

    if (Pstream::parRun())
    {
        dirDetailBlade = mesh_.time().path()/"../postProcessing"
            / name_ / mesh_.time().timeName();
    }
    else
    {
        dirDetailBlade = mesh_.time().path()/"postProcessing"
            / name_ / mesh_.time().timeName();
    }

    if (not isDir(dirDetailBlade))
    {
        mkDir(dirDetailBlade);
    }


    if (Pstream::master())
    {
        Pout << "createOutputFile() blade detail, processor: " << currentProcsID_ << endl;
        Pout << "createOutputFile() blade detail, bladeNumberofActuatorPoints: " 
            << bladeActuatorPointsPositions_.size() << endl;

        outputFileDetailBlade_ = new OFstream(dirDetailBlade/"detailAboutBlade.csv");
        *outputFileDetailBlade_<< "Time" << endl;
        *outputFileDetailBlade_<< "blade number, element number, processorID, (x, y, z), (avg_u, avg_v, avg_w), (u, v, w)," 
            << "(fx, fy, fz), bladeLocalCl, bladeLocalCd, bladePitch:" << endl;
    }


    // rotor detail

    fileName dirDetailRotor;

    if (Pstream::parRun())
    {
        dirDetailRotor = mesh_.time().path()/"../postProcessing"
            / name_ / mesh_.time().timeName();
    }
    else
    {
        dirDetailRotor = mesh_.time().path()/"postProcessing"
            / name_ / mesh_.time().timeName();
    }

    if (not isDir(dirDetailRotor))
    {
        mkDir(dirDetailRotor);
    }


    if (Pstream::master())
    {
        Pout << "createOutputFile() rotor detail, processor: " << currentProcsID_ << endl;
        Pout << "createOutputFile() rotor detail, rotorNumberofActuatorPoints: " 
            << rotorActuatorPointsPositions_.size() << endl;

        outputFileDetailRotor_ = new OFstream(dirDetailRotor/"detailAboutRotor.csv");
        *outputFileDetailRotor_<< "Time" << endl;
        *outputFileDetailRotor_<< "element number, processorID, (x, y, z), (u, v, w), (fx, fy, fz):" << endl;
    }


    // General

    fileName dirGeneral;

    if (Pstream::parRun())
    {
        dirGeneral = mesh_.time().path()/"../postProcessing"
            / name_ / mesh_.time().timeName();
    }
    else
    {
        dirGeneral = mesh_.time().path()/"postProcessing"
            / name_ / mesh_.time().timeName();
    }

    if (not isDir(dirGeneral))
    {
        mkDir(dirGeneral);
    }


    if (Pstream::master())
    {
        Pout << "createOutputFile() general, processor: " << currentProcsID_ << endl;
        Pout << "createOutputFile() general, rotorNumberofActuatorPoints: " << rotorForceExertingActiveNum_ << endl;

        outputFileGeneral_ = new OFstream(dirGeneral/"generalAboutMRSL.csv");
        *outputFileGeneral_<< "Time and then rotorCtStar, (Ux, Uy, Uz)_DiskAVG, (fx, fy, fz)_rotor, (fx, fy, fz)_blade, (fx, fy, fz)_rotor_calc, (fx, fy, fz)_blade_calc, power_rotor_calc, power_blade_calc:" << endl;
    }

    
}



Foam::scalar Foam::fv::flyingActuationDiskSource::interpolate
(
    scalar xNew,
    scalarList& xOld,
    scalarList& yOld
)
{
    label index = 0;
    label indexP = 0;
    label indexM = 0;
    scalar error = 1.0E30;
    forAll(xOld, i)
    {
        scalar diff = mag(xNew - xOld[i]);
        if (diff < error)
        {
            index = i;
            error = diff;
        }
    }
    if (xNew < xOld[index])
    {
        if (index == 0)
        {
            indexP = 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index;
            indexM = indexP - 1;
        }
        return yOld[indexM]
               + ((yOld[indexP]
               - yOld[indexM])/(xOld[indexP]
               - xOld[indexM]))*(xNew - xOld[indexM]);
    }
    else if (xNew > xOld[index])
    {
        if (index == xOld.size() - 1)
        {
            indexP = xOld.size() - 1;
            indexM = indexP - 1;
        }
        else
        {
            indexP = index + 1;
            indexM = indexP - 1;
        }
        return yOld[indexM] + ((yOld[indexP]
               - yOld[indexM])/(xOld[indexP]
               - xOld[indexM]))*(xNew - xOld[indexM]);
    }
    else if (xNew == xOld[index])
    {
        return yOld[index];
    }
    else
    {
        return 0.0;
    }
}

// ************************************************************************* //
