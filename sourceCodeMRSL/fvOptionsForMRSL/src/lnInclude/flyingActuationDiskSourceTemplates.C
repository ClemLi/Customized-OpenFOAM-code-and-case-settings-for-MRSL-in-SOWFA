/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------

\*---------------------------------------------------------------------------*/

#include "flyingActuationDiskSource.H"

// * * * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * //


template<class AlphaFieldType, class RhoFieldType>
void Foam::fv::flyingActuationDiskSource::makeWindTurbineWakesFly
(
    const AlphaFieldType& alpha,
    const RhoFieldType& rho,
    fvMatrix<vector>& eqn,
    interpolationCellPoint<vector> UInterp
)
{
    /////////////////////////////////// General initialization ///////////////////////////////////

    label sampleCellI = -1;


    // Zero out the forces
    localForceRotor_ = vector::zero;
    localForceRotorCalc_ = vector::zero;
    localForceBlade_ = vector::zero;
    localForceBladeCalc_ = vector::zero;

    rotorPowerCalc_ = 0.0;
    bladePowerCalc_ = 0.0;

    rotorALForceList_ = vector::zero;
    bladeALForceList_ = vector::zero;

    // Caution: values of forceField_ are only valid when the time step count is an integer multiples of writeInterval_
    if ( ( rotorForceRecievingActiveNum_ + bladeForceRecievingActiveNum_ ) > 0 && (writeIntervalTrack_+1) == writeInterval_ )
    {
        // Zero out force vector and field
        forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);

        // Check dimensions of force field and correct if necessary
        if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
        {
            forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
        }
    }
    else if ( currentTimeStep_ < 2.0 )
    {
        // Zero out force vector and field
        forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);

        // Check dimensions of force field and correct if necessary
        if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
        {
            forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
        }
    }


    const vectorField& U = eqn.psi();

    scalar currentTime = mesh_.time().value();


    /////////////////////////////////// Pitching the blades to the set AOA ///////////////////////////////////


    if( timeStartAdjustPitch_ > 0 &&  currentTime >= timeStartAdjustPitch_ )
    {

        for( label bladeIndex = 0; bladeIndex < bladeNumber_; bladeIndex++ )
        {
            scalar misTargetAOANow = bladeAOAMidofThatBlade_[bladeIndex] - bladeAOATarget_[bladeIndex];

            misTargetAOAFiltered_[bladeIndex] = 
                (misTargetAOANow + misTargetAOAFiltered_[bladeIndex]) / bladePitchTimeScale_ 
                * exp( - 1 / bladePitchTimeScale_ *  mesh_.time().deltaT().value() );

            scalar deltaPitch = bladePitchSpeed_ * mesh_.time().deltaT().value();

            if( misTargetAOAFiltered_[bladeIndex] > ( bladePitchTol_ * deltaPitch ) )
            {
            for (label indexofAL = 0; indexofAL < bladeActuatorPointsPositions_.size(); indexofAL++ )
                {
                    if(bladeNumIndicator_[indexofAL] == bladeIndex)
                    {
                        bladePitch_[indexofAL] = bladePitch_[indexofAL] - deltaPitch;
                    }
                }
            }

            if(misTargetAOAFiltered_[bladeIndex] < ( -1 * bladePitchTol_ * deltaPitch ) )
            {
                for (label indexofAL = 0; indexofAL < bladeActuatorPointsPositions_.size(); indexofAL++ )
                {
                     if(bladeNumIndicator_[indexofAL] == bladeIndex)
                    {
                        bladePitch_[indexofAL] = bladePitch_[indexofAL] + deltaPitch;
                    }
                }
            }
        }
   }


    /////////////////////////////////// Initialization for the rotor ///////////////////////////////////

    label rotorOrBlade = 0;  // 0 stands for rotor, 1 stands for blade

    inflowVelocityRotorList_ = vector::zero;
    vector u_disk = Foam::Zero;


    /////////////////////////////////// Calculating the forces on the rotor-AL points ///////////////////////////////////


    if ( rotorForceExertingActiveNum_ > 0 )
    {

        if ( currentTimeStep_ == timeStepToCheck_ )
        {
            Pout << name_ << " Exerting: rotorForceExertingActiveNum_: " << rotorForceExertingActiveNum_ << endl;
            Pout << endl;

            rotorForceActiveNumChecker_ = 0;
        }

        for 
        (
            label actuatorPointIndex = 0; 
            actuatorPointIndex < rotorActuatorPointsPositions_.size();
            actuatorPointIndex++
        )
        {

            if ( rotorExertingActiveList_[actuatorPointIndex] == 1 )
            {

                if ( currentTimeStep_ == timeStepToCheck_ )
                {
                    Pout << "Rotor: before getting the force at point " << actuatorPointIndex
                        << " for " << name_ << ", which is: "
                        << rotorActuatorPointsPositions_[actuatorPointIndex] << endl;

                    rotorForceActiveNumChecker_++;
                }

                vector cellCenterNow = rotorActuatorPointsPositions_[actuatorPointIndex];

                vector rotorForceVectorToApply = Foam::Zero;

                sampleCellI = rotorCellID_[actuatorPointIndex];
                vector inflowVelocityRotor = UInterp.interpolate( cellCenterNow, sampleCellI );

                inflowVelocityRotorList_[actuatorPointIndex] = inflowVelocityRotor;

                u_disk = u_disk + inflowVelocityRotor / rotorActuatorPointsPositions_.size();

                rotorForceVectorToApply = 0.5 * (inflowVelocityRotor &  diskDir_) * (inflowVelocityRotor & diskDir_)
                    * rotorFrontalAreaPerCell_[actuatorPointIndex] * rotorCtStar_ * diskDir_;

                localForceRotor_ = localForceRotor_ + rotorForceVectorToApply;

                rotorALForceList_[actuatorPointIndex] = rotorForceVectorToApply;
            }
        }

        if (currentTimeStep_ == timeStepToCheck_)
        {
            if (rotorForceActiveNumChecker_ != rotorForceExertingActiveNum_ )
            {
                FatalErrorIn("rotorForceActiveNumChecker_ eror")
                    << "rotorForceActiveNumChecker_ for " << name_
                    << " in this processor is " << rotorForceActiveNumChecker_ << ", but rotorForceExertingActiveNum_ "
                    << " is " << rotorForceExertingActiveNum_ << endl
                    << "note that we find that variables in OpenFOAM sometimes tend to be volatile"
                    << abort(FatalError);
            }
            else
            {
                Pout << "rotorForceActiveNumChecker_ for " << name_
                    << " and rotorForceExertingActiveNum_ in this processor are both " << rotorForceActiveNumChecker_
                    << endl;
            }
        }
    }


    /////////////////////////////////// Broadcast the rotor-AL forces and rotor-related parameters ///////////////////////////////////

    // sum onto master
    Pstream::listCombineGather(rotorALForceList_, plusEqOp<vector>());
    Pstream::listCombineGather(inflowVelocityRotorList_, plusEqOp<vector>());

    // broadcast summed lists to everyone
    // inflowVelocityRotorList_ is only used by the master processor, no need to be broadcasted
    Pstream::listCombineScatter(rotorALForceList_);

    /////////////////////////////////// Recieving the smeared body force fields of rotor ///////////////////////////////////

    if ( rotorForceRecievingActiveNum_ > 0 )
    {
        if ( currentTimeStep_ == timeStepToCheck_ )
        {
            Pout << "Recieve: rotorForceRecievingActiveNum_:" << " for " << name_ 
                << ": " << rotorForceRecievingActiveNum_ << endl;

            rotorForceActiveNumChecker_ = 0;
        }

        for 
        (
            label actuatorPointIndex = 0; 
            actuatorPointIndex < rotorActuatorPointsPositions_.size();
            actuatorPointIndex++
        )
        {

            if ( rotorLabelListtoSmear_[actuatorPointIndex][0] > -1 )
            {
                if ( currentTimeStep_ == timeStepToCheck_ )
                {
                    Pout << name_ << " Rotor: before applying applyForceField(), " 
                        << "rotorLabelListtoSmear_[" << actuatorPointIndex << "][0]: " 
                        << rotorLabelListtoSmear_[actuatorPointIndex][0] << endl;
                    Pout << "rotorLabelListtoSmear_[" << actuatorPointIndex << "].size(): " 
                        << rotorLabelListtoSmear_[actuatorPointIndex].size() << endl;

                    rotorForceActiveNumChecker_++;
                }

                vector cellCenterNow = rotorActuatorPointsPositions_[actuatorPointIndex];

                vector rotorForceVectorToApply = rotorALForceList_[actuatorPointIndex];

                applyForceField(eqn, U, cellCenterNow, rotorForceVectorToApply, 
                    rotorLabelListtoSmear_[actuatorPointIndex], rotorSmearFactorListList_[actuatorPointIndex],
                    rotorOrBlade);
            }

        }

        if (currentTimeStep_ == timeStepToCheck_)
        {
            if (rotorForceActiveNumChecker_ != rotorForceRecievingActiveNum_ )
            {
                FatalErrorIn("rotorForceActiveNumChecker_ eror")
                    << "rotorForceActiveNumChecker_ for " << name_
                    << " in this processor is " << rotorForceActiveNumChecker_ << ", but rotorForceRecievingActiveNum_ "
                    << " is " << rotorForceRecievingActiveNum_ << endl
                    << "note that we find that variables in OpenFOAM sometimes tend to be volatile"
                    << abort(FatalError);
            }
            else
            {
                Pout << "rotorForceActiveNumChecker_ for " << name_
                    << " and rotorForceRecievingActiveNum_ in this processor are both " << rotorForceActiveNumChecker_
                    << endl;
            }
        }

    }


    /////////////////////////////////// Initialization for the blades ///////////////////////////////////


    rotorOrBlade = 1;

    //bladeFlowAngle_ = scalarList(bladeActuatorPointsPositions_.size(), 0.0);
    bladeFlowAngle_ = 0.0;  // still a scalarList
    bladeAOA_ = 0.0;

    bladeClList_ = 0.0;
    bladeCdList_ = 0.0;

    inflowVelocityList_ = vector::zero;
    inflowVelocityOnPointList_ = vector::zero;


    /////////////////////////////////// Get the inflow velocity for the blade-AL points ///////////////////////////////////


    for 
    (
        label actuatorPointIndex = 0; 
        actuatorPointIndex < bladeActuatorPointsPositions_.size();
        actuatorPointIndex++
    )
    {
        vector cellCenterNow = bladeActuatorPointsPositions_[actuatorPointIndex];

        scalar velSamplingRadius = bladeSmearingLength_[actuatorPointIndex] * velSamplingRadiusCoeff_;

        inflowVelocity_ = vector::zero;

        // Ring-line sampling
        if( velSamplingRadius > 0.0 ) 
        {
            if ( ringVelocityGlobalChecker_ > 0.0 )
            {
                for( label velSampleIndex = 0; velSampleIndex < velSamplingNumber_; velSampleIndex++ )
                {
                    // distribute the points evenly in terms of angular distance
                    scalar pointAngle = Foam::constant::mathematical::pi * 2.0 * velSampleIndex / velSamplingNumber_;
                    scalar freeStreamDist = velSamplingRadius * Foam::cos(pointAngle);
                    scalar verticalDist = velSamplingRadius * Foam::sin(pointAngle);
                    vector samplePoint = bladeActuatorPointsPositions_[actuatorPointIndex]
                        + freeStreamDist * freeStreamDir_ + verticalDist * (freeStreamDir_ ^ lateralDirection_);

                    // The cell to look into has already been pre-sampled, if the point lies in other processor, it will be -1
                    sampleCellI = inflowVelocityRingSampleList_[actuatorPointIndex][velSampleIndex];

                    if ( sampleCellI > -1 )
                    {
                        if ( currentTimeStep_ == timeStepToCheck_ )
                        {
                            Pout << "checking: inflowVelocityRingSampleList_[" << actuatorPointIndex << "][" << velSampleIndex << "]: " 
                                << inflowVelocityRingSampleList_[actuatorPointIndex][velSampleIndex] << endl;
                        }

                        inflowVelocity_ = inflowVelocity_ + UInterp.interpolate( samplePoint, sampleCellI ) / velSamplingNumber_;
                    }
                }

            }

        }
        
        // Single point sampling

        inflowVelocityOnPoint_ = vector::zero;

        sampleCellI = bladeCellID_[actuatorPointIndex];

        if ( sampleCellI > -1 )
        {
            inflowVelocityOnPoint_ = inflowVelocityOnPoint_ + UInterp.interpolate( cellCenterNow, sampleCellI );
            if ( currentTimeStep_ == timeStepToCheck_ )
            {
                Pout << "checking: sampleCellI: " << sampleCellI << endl;
                Pout << "checking: inflowVelocityOnPoint_: " << inflowVelocityOnPoint_ << endl;
            }
        }

        if ( velSamplingRadius <= 0.0 )
        {
            inflowVelocity_ = inflowVelocityOnPoint_;
        } 

        inflowVelocityList_[actuatorPointIndex] = inflowVelocity_;
        inflowVelocityOnPointList_[actuatorPointIndex] = inflowVelocityOnPoint_;
    }


    // scatter inflowVelocities

    // sum onto master
    Pstream::listCombineGather(inflowVelocityList_, plusEqOp<vector>());
    Pstream::listCombineGather(inflowVelocityOnPointList_, plusEqOp<vector>());

    // broadcast summed lists to everyone
    // inflowVelocityOnPointList_ is only used by the master, no need to be broadcasted
    Pstream::listCombineScatter(inflowVelocityList_);


    /////////////////////////////////// Calculating the forces on the blades-AL points ///////////////////////////////////
 
    
    if ( bladeForceExertingActiveNum_ > 0 )
    {

        if ( currentTimeStep_ == timeStepToCheck_ )
        {
            Pout << name_ << " Exerting: bladeForceExertingActiveNum_" << " for " << name_ 
                << ": " << bladeForceExertingActiveNum_ << endl;

            bladeForceActiveNumChecker_ = 0;
        }

        for 
        (
            label actuatorPointIndex = 0; 
            actuatorPointIndex < bladeActuatorPointsPositions_.size();
            actuatorPointIndex++
        )
        {

            if ( bladeExertingActiveList_[actuatorPointIndex] == 1 )
            {

                if ( currentTimeStep_ == timeStepToCheck_ )
                {
                    Pout << "Blade: before getting the force at point " << actuatorPointIndex
                        << " for " << name_ << ", which is " << bladeActuatorPointsPositions_[actuatorPointIndex] << endl;

                    bladeForceActiveNumChecker_++;
                }

                inflowVelocity_ = inflowVelocityList_[actuatorPointIndex];

                scalar samplingFlowAngle;
                samplingFlowAngle = Foam::atan(inflowVelocity_.z() / (inflowVelocity_ & freeStreamDir_) );
                bladeFlowAngle_[actuatorPointIndex] = samplingFlowAngle;

                bladeAOA_[actuatorPointIndex] = 
                    ( samplingFlowAngle * 180 / Foam::constant::mathematical::pi + 
                    bladeTwist_[actuatorPointIndex] + bladePitch_[actuatorPointIndex] );

                scalar bladeLocalCl = interpolate(bladeAOA_[actuatorPointIndex], alphaList_, clList_);

                scalar bladeLocalCd = interpolate(bladeAOA_[actuatorPointIndex], alphaList_, cdList_);

                bladeClList_[actuatorPointIndex] = bladeLocalCl;
                bladeCdList_[actuatorPointIndex] = bladeLocalCd;


                vector bladeLiftVectorToApply = 0.5 *
                    ( (inflowVelocity_ & freeStreamDir_) * (inflowVelocity_ & freeStreamDir_) + inflowVelocity_.z() * inflowVelocity_.z() ) *
                    bladeLocalCl * bladeChord_ * bladeSpanPerAL_ * 
                    ( (inflowVelocity_/ mag(inflowVelocity_) ) ^ lateralDirection_ );

                vector bladeDragVectorToApply = 0.5 * 
                    ( (inflowVelocity_ & freeStreamDir_) * (inflowVelocity_ & freeStreamDir_) + inflowVelocity_.z() * inflowVelocity_.z() ) *
                    bladeLocalCd * bladeChord_ * bladeSpanPerAL_ * 
                    ( inflowVelocity_ / mag(inflowVelocity_) );

                bladeDragVectorToApply.y() = 0.0;   // To mute spawise directional force

                vector bladeForceVectorToApply = bladeLiftVectorToApply + bladeDragVectorToApply;

                localForceBlade_ = localForceBlade_ + bladeForceVectorToApply;

                bladeALForceList_[actuatorPointIndex] = bladeForceVectorToApply;
            }
        }


        if (currentTimeStep_ == timeStepToCheck_)
        {
            if (bladeForceActiveNumChecker_ != bladeForceExertingActiveNum_ )
            {
                FatalErrorIn("bladeForceActiveNumChecker_ eror")
                    << "bladeForceActiveNumChecker_ for " << name_
                    << " in this processor is " << bladeForceActiveNumChecker_ << ", but bladeForceExertingActiveNum_ "
                    << " is " << bladeForceExertingActiveNum_ << endl
                    << "note that we find that variables in OpenFOAM sometimes tend to be volatile"
                    << abort(FatalError);
            }
            else
            {
                Pout << "bladeForceActiveNumChecker_ for " << name_
                    << " and bladeForceExertingActiveNum_ are in this processor both " << bladeForceActiveNumChecker_
                    << endl;
            }
        }
    }

    /////////////////////////////////// Broadcast the blade-AL forces and blade-related parameters ///////////////////////////////////


    // sum onto master
    Pstream::listCombineGather(bladeALForceList_, plusEqOp<vector>());
    Pstream::listCombineGather(bladeClList_, plusEqOp<scalar>());
    Pstream::listCombineGather(bladeCdList_, plusEqOp<scalar>());
    Pstream::listCombineGather(bladeAOA_, plusEqOp<scalar>());

    // broadcast summed lists to everyone
    // Cl and Cd are only used by the master processor, no need to be broadcasted
    Pstream::listCombineScatter(bladeALForceList_);
    Pstream::listCombineScatter(bladeAOA_);


    // To get the AOA at mid span of each blade
    for(label bladeIndex = 0; bladeIndex < bladeNumber_; bladeIndex++)
    {
        bladeAOAMidofThatBlade_[bladeIndex] = bladeAOA_[indexForBladeMid_[bladeIndex]];
    }


    /////////////////////////////////// Recieving the smeared body force fields of blades ///////////////////////////////////

    if ( bladeForceRecievingActiveNum_ > 0 )
    {

        if ( currentTimeStep_ == timeStepToCheck_ )
        {
            Pout << "Recieve: bladeForceRecievingActiveNum_: " << bladeForceRecievingActiveNum_ << endl;
            Pout << endl;

            bladeForceActiveNumChecker_ = 0;
        }

        for 
        (
            label actuatorPointIndex = 0; 
            actuatorPointIndex < bladeActuatorPointsPositions_.size();
            actuatorPointIndex++
        )
        {
            if ( bladeLabelListtoSmear_[actuatorPointIndex][0] > -1 )
            {

                if ( currentTimeStep_ == timeStepToCheck_ )
                {
                    Pout << name_ << " Blade: before applying applyForceField() at point, " 
                        << "bladeLabelListtoSmear_[" << actuatorPointIndex << "][0]: " 
                        << bladeLabelListtoSmear_[actuatorPointIndex][0] << endl;
                    Pout << "bladeLabelListtoSmear_[" << actuatorPointIndex << "].size(): " 
                        << bladeLabelListtoSmear_[actuatorPointIndex].size() << endl;

                    bladeForceActiveNumChecker_++;
                }

                vector cellCenterNow = bladeActuatorPointsPositions_[actuatorPointIndex];

                vector bladeForceVectorToApply = bladeALForceList_[actuatorPointIndex];

                applyForceField(eqn, U, cellCenterNow, bladeForceVectorToApply, 
                    bladeLabelListtoSmear_[actuatorPointIndex], bladeSmearFactorListList_[actuatorPointIndex],
                    rotorOrBlade);
            }    
        }

        if (currentTimeStep_ == timeStepToCheck_)
        {
            if (bladeForceActiveNumChecker_ != bladeForceRecievingActiveNum_ )
            {
                FatalErrorIn("bladeForceActiveNumChecker_ eror")
                    << "bladeForceActiveNumChecker_ for " << name_
                    << " in this processor is " << bladeForceActiveNumChecker_ << ", but bladeForceRecievingActiveNum_ "
                    << " is " << bladeForceRecievingActiveNum_ << endl
                    << "note that we find that variables in OpenFOAM sometimes tend to be volatile"
                    << abort(FatalError);
            }
            else
            {
                Pout << "bladeForceActiveNumChecker_ for " << name_
                    << " and bladeForceRecievingActiveNum_ are in this processor both " << bladeForceActiveNumChecker_
                    << endl;
            }
        }

    }

    /////////////////////////////////// Collect and add up the forces in all processors ///////////////////////////////////

    outputFileTime_[writeIntervalTrack_] = mesh_.time().value();


    outputFileDataGeneral_[writeIntervalTrack_][0]  = localForceRotor_.x();
    outputFileDataGeneral_[writeIntervalTrack_][1]  = localForceRotor_.y();
    outputFileDataGeneral_[writeIntervalTrack_][2]  = localForceRotor_.z();

    outputFileDataGeneral_[writeIntervalTrack_][3]  = localForceRotorCalc_.x();
    outputFileDataGeneral_[writeIntervalTrack_][4]  = localForceRotorCalc_.y();
    outputFileDataGeneral_[writeIntervalTrack_][5]  = localForceRotorCalc_.z();


    outputFileDataGeneral_[writeIntervalTrack_][6]  = localForceBlade_.x();
    outputFileDataGeneral_[writeIntervalTrack_][7]  = localForceBlade_.y();
    outputFileDataGeneral_[writeIntervalTrack_][8]  = localForceBlade_.z();

    outputFileDataGeneral_[writeIntervalTrack_][9]  = localForceBladeCalc_.x();
    outputFileDataGeneral_[writeIntervalTrack_][10] = localForceBladeCalc_.y();
    outputFileDataGeneral_[writeIntervalTrack_][11] = localForceBladeCalc_.z();


    outputFileDataGeneral_[writeIntervalTrack_][12] = rotorPowerCalc_;

    outputFileDataGeneral_[writeIntervalTrack_][13] = bladePowerCalc_;


    outputFileDataGeneral_[writeIntervalTrack_][14] = u_disk.x();
    outputFileDataGeneral_[writeIntervalTrack_][15] = u_disk.y();
    outputFileDataGeneral_[writeIntervalTrack_][16] = u_disk.z();


    // sum onto master
    Pstream::listCombineGather(outputFileDataGeneral_[writeIntervalTrack_], plusEqOp<scalar>());


    /////////////////////////////////// Log the rotor-AL information to buffer for output ///////////////////////////////////

    if ( Pstream::master() )
    {

        for 
        (
            label actuatorPointIndex = 0; 
            actuatorPointIndex < rotorActuatorPointsPositions_.size();
            actuatorPointIndex++
        )
        {
            outputFileRotorInflowVelX_[writeIntervalTrack_][actuatorPointIndex]  = inflowVelocityRotorList_[actuatorPointIndex].x();
            outputFileRotorInflowVelY_[writeIntervalTrack_][actuatorPointIndex]  = inflowVelocityRotorList_[actuatorPointIndex].y();
            outputFileRotorInflowVelZ_[writeIntervalTrack_][actuatorPointIndex]  = inflowVelocityRotorList_[actuatorPointIndex].z();

            outputFileRotorALForceX_[writeIntervalTrack_][actuatorPointIndex]  = rotorALForceList_[actuatorPointIndex].x();
            outputFileRotorALForceY_[writeIntervalTrack_][actuatorPointIndex]  = rotorALForceList_[actuatorPointIndex].y();
            outputFileRotorALForceZ_[writeIntervalTrack_][actuatorPointIndex]  = rotorALForceList_[actuatorPointIndex].z();

        }
    }


    /////////////////////////////////// Log the rotor-AL information to buffer for output ///////////////////////////////////

    if ( Pstream::master() )
    {
        for 
        (
            label actuatorPointIndex = 0; 
            actuatorPointIndex < bladeActuatorPointsPositions_.size();
            actuatorPointIndex++
        )
        {

            outputFileBladeInflowVelX_[writeIntervalTrack_][actuatorPointIndex] = inflowVelocityList_[actuatorPointIndex].x();
            outputFileBladeInflowVelY_[writeIntervalTrack_][actuatorPointIndex] = inflowVelocityList_[actuatorPointIndex].y();
            outputFileBladeInflowVelZ_[writeIntervalTrack_][actuatorPointIndex] = inflowVelocityList_[actuatorPointIndex].z();

            outputFileBladeInflowVelOnPointX_[writeIntervalTrack_][actuatorPointIndex] = inflowVelocityOnPointList_[actuatorPointIndex].x();
            outputFileBladeInflowVelOnPointY_[writeIntervalTrack_][actuatorPointIndex] = inflowVelocityOnPointList_[actuatorPointIndex].y();
            outputFileBladeInflowVelOnPointZ_[writeIntervalTrack_][actuatorPointIndex] = inflowVelocityOnPointList_[actuatorPointIndex].z();

            outputFileBladeALForceX_[writeIntervalTrack_][actuatorPointIndex] = bladeALForceList_[actuatorPointIndex].x();
            outputFileBladeALForceY_[writeIntervalTrack_][actuatorPointIndex] = bladeALForceList_[actuatorPointIndex].y();
            outputFileBladeALForceZ_[writeIntervalTrack_][actuatorPointIndex] = bladeALForceList_[actuatorPointIndex].z();

            outputFileBladeCl_[writeIntervalTrack_][actuatorPointIndex] = bladeClList_[actuatorPointIndex];
            outputFileBladeCd_[writeIntervalTrack_][actuatorPointIndex] = bladeCdList_[actuatorPointIndex];
            outputFileBladePitch_[writeIntervalTrack_][actuatorPointIndex] = bladePitch_[actuatorPointIndex];

        }
    }

    /////////////////////////////////// Write the out put files ///////////////////////////////////


    if ( Pstream::master() && outPutTimeStart_ <= currentTime && (writeIntervalTrack_+1) == writeInterval_ )
    {

        //// rotor
        for(label writeIdx = 0; writeIdx < writeInterval_; writeIdx++)
        {

            *outputFileDetailRotor_ << outputFileTime_[writeIdx] << "," << '\n'; 

            for 
            (
                label actuatorPointIndex = 0; 
                actuatorPointIndex < rotorActuatorPointsPositions_.size();
                actuatorPointIndex++
            )
            {
                *outputFileDetailRotor_ << actuatorPointIndex << ","

                << round( rotorALWhichProcs_[actuatorPointIndex] ) << ","

                << "(" << rotorActuatorPointsPositions_[actuatorPointIndex].x() << "," 
                    << rotorActuatorPointsPositions_[actuatorPointIndex].y() << "," 
                    << rotorActuatorPointsPositions_[actuatorPointIndex].z() << ")" << ","

                << "(" << outputFileRotorInflowVelX_[writeIdx][actuatorPointIndex] << "," 
                    << outputFileRotorInflowVelY_[writeIdx][actuatorPointIndex] << "," 
                    << outputFileRotorInflowVelZ_[writeIdx][actuatorPointIndex] << ")" << ","

                << "(" << outputFileRotorALForceX_[writeIdx][actuatorPointIndex] * -1.0 << "," 
                    << outputFileRotorALForceY_[writeIdx][actuatorPointIndex] * -1.0 << "," 
                    << outputFileRotorALForceZ_[writeIdx][actuatorPointIndex] * -1.0 << ")" << ","

                << '\n';
            }

            *outputFileDetailRotor_ << '\n';
        }


        //// blade

        for(label writeIdx = 0; writeIdx < writeInterval_; writeIdx++)
        {

            *outputFileDetailBlade_ <<outputFileTime_[writeIdx] << "," << '\n';

            for 
            (
                label actuatorPointIndex = 0; 
                actuatorPointIndex < bladeActuatorPointsPositions_.size();
                actuatorPointIndex++
            )
            {
                *outputFileDetailBlade_ << bladeNumIndicator_[actuatorPointIndex] << ","

                << actuatorPointIndex << ","

                << round( bladeALWhichProcs_[actuatorPointIndex] ) << ","

                << "(" << bladeActuatorPointsPositions_[actuatorPointIndex].x() << "," 
                    << bladeActuatorPointsPositions_[actuatorPointIndex].y() << "," 
                    << bladeActuatorPointsPositions_[actuatorPointIndex].z() << ")" << ","

                << "(" << outputFileBladeInflowVelX_[writeIdx][actuatorPointIndex] << "," 
                    << outputFileBladeInflowVelY_[writeIdx][actuatorPointIndex] << "," 
                    << outputFileBladeInflowVelZ_[writeIdx][actuatorPointIndex] << ")" << ","

                << "(" << outputFileBladeInflowVelOnPointX_[writeIdx][actuatorPointIndex] << "," 
                    << outputFileBladeInflowVelOnPointY_[writeIdx][actuatorPointIndex] << "," 
                    << outputFileBladeInflowVelOnPointZ_[writeIdx][actuatorPointIndex] << ")" << ","

                << "(" << outputFileBladeALForceX_[writeIdx][actuatorPointIndex] * -1.0 << "," 
                    << outputFileBladeALForceY_[writeIdx][actuatorPointIndex] * -1.0 << "," 
                    << outputFileBladeALForceZ_[writeIdx][actuatorPointIndex] * -1.0 << ")" << ","

                << outputFileBladeCl_[writeIdx][actuatorPointIndex] << ","

                << outputFileBladeCd_[writeIdx][actuatorPointIndex] << ","

                << outputFileBladePitch_[writeIdx][actuatorPointIndex] << ","

                << '\n';

            }

            *outputFileDetailBlade_ << '\n';
        }


        //// general

        for(label writeIdx = 0; writeIdx < writeInterval_; writeIdx++)
        {

            *outputFileGeneral_ << outputFileTime_[writeIdx] << ","

            << rotorCtStar_ << ","

            // u_disk
            << "(" << outputFileDataGeneral_[writeIdx][14] << "," 
                << outputFileDataGeneral_[writeIdx][15] << "," << outputFileDataGeneral_[writeIdx][16] << ")" << ","

            // localForceRotor_
            << "(" << outputFileDataGeneral_[writeIdx][0] << "," 
                << outputFileDataGeneral_[writeIdx][1] << "," << outputFileDataGeneral_[writeIdx][2] << ")" << ","

            // localForceBlade_
            << "(" << outputFileDataGeneral_[writeIdx][6] << "," 
                << outputFileDataGeneral_[writeIdx][7] << "," << outputFileDataGeneral_[writeIdx][8] << ")" << ","

            // localForceRotorCalc_
            << "(" << outputFileDataGeneral_[writeIdx][3] << "," 
                << outputFileDataGeneral_[writeIdx][4] << "," << outputFileDataGeneral_[writeIdx][5] << ")" << ","

            // localForceBladeCalc_
            << "(" << outputFileDataGeneral_[writeIdx][9] << "," 
                << outputFileDataGeneral_[writeIdx][10] << "," << outputFileDataGeneral_[writeIdx][11] << ")" << ","

            // rotorPowerCalc_
            << outputFileDataGeneral_[writeIdx][12] << ","

            // bladePowerCalc_
            << outputFileDataGeneral_[writeIdx][13] << ","

            << '\n';

        }

    }


    // update writeIntervalTrack_


    if( (writeIntervalTrack_+1) == writeInterval_ )
    {

        if( Pstream::master() )
        {
            Pout << endl << name_ << ": Flushing the outout files ..." << endl << endl;

            outputFileDetailRotor_->flush();  // flush once per dump
            outputFileDetailBlade_->flush();  // flush once per dump
            outputFileGeneral_->flush();  // flush once per dump
        }

        writeIntervalTrack_ = 0;

    }
    else
    {
        writeIntervalTrack_ += 1;
    }


}



void Foam::fv::flyingActuationDiskSource::applyForceField
(
    fvMatrix<vector>& eqn,
    const vectorField& U,
    vector cellCenterNow,
    vector forceVectorToApply,
    labelList domainCells,  //Cells selected to smear
    scalarList smearFactorForDomainCells,
    label rotorOrBlade
)
{

    // Apply force to the cells within the element's sphere of influence
    for (label elementNow = 0; elementNow < domainCells.size(); elementNow++)
    {
        label cellI = domainCells[elementNow];

        scalar factor = smearFactorForDomainCells[elementNow];

        vector forceSmearedApplyToCell = forceVectorToApply*factor*mesh_.V()[cellI]; // become force per density after smearing
        
        // 0 stands for rotor, 1 stands for blade
        if(rotorOrBlade == 0)
        {
            localForceRotorCalc_ -= forceSmearedApplyToCell;
            rotorPowerCalc_ += forceSmearedApplyToCell & U[cellI];   // The work actuator posed on the the flow
        }

        // 0 stands for rotor, 1 stands for blade
        if(rotorOrBlade == 1)
        {
            localForceBladeCalc_ -= forceSmearedApplyToCell;
            bladePowerCalc_ += forceSmearedApplyToCell & U[cellI];
        }

        forceField_[cellI] -= forceVectorToApply*factor;  //force per volume per density
        eqn.source()[cellI] += forceSmearedApplyToCell;    // Informing the NSE with the momentum source

    }
}


// ************************************************************************* //
