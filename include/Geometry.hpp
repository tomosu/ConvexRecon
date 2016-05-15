#pragma once

#include "RealTypes.hpp"
#include "libini.h"
#include <cassert>
#include <string>
#include <iostream>

#define MBF 256

struct Geometry {

  enum class BeamType{
    PARALLEL_BEAM,
      CONE_BEAM,
  };

  enum class DetectorType{
    HORIZONTAL,
  };

  //general
  Real_t sid; //SourceImageDistance
  Real_t sod; //SourceObjectDistance
  BeamType beamType;
  Real2_t gantryCenter;

  //source
  Real2_t sourceOrigin;

  //detector(projection)
  DetectorType detectorType;
  Real2_t detectorCenterPos;
  Real_t detectorCenterPixel;
  int detectorPixelNum;
  Real_t detectorPitch;
  int projectionCount;

  //volume(reconstruction)
  Int2_t volumeSize;
  Real2_t volumeCenter;
  Real2_t volumePitch;


  Geometry(char *filepath){

    INICFG ini;

    //general
    ini.AddEntry("SID",sid); //SourceImageDistance
    ini.AddEntry("SOD",sod); //SourceObjectDistance
    char btype[MBF] ="BEAM";
    ini.AddEntry("BEAM_TYPE", btype, MBF);
    ini.AddEntry("GANTRY_CENTER_X",gantryCenter.x);
    ini.AddEntry("GANTRY_CENTER_Y",gantryCenter.y);

    //detector(projection)
    char dtype[MBF] ="DETECTOR";
    ini.AddEntry("DETECTOR_TYPE",dtype, MBF);

    ini.AddEntry("DETECTOR_CENTER_PIXEL",detectorCenterPixel);
    ini.AddEntry("DETECTOR_PIXEL_NUM",detectorPixelNum);
    ini.AddEntry("DETECTOR_PITCH",detectorPitch);
    ini.AddEntry("PROJECTION_COUNT",projectionCount);

    //volume(reconstruction)
    ini.AddEntry("VOLUME_SIZE_X",volumeSize.x);
    ini.AddEntry("VOLUME_SIZE_Y",volumeSize.y);
    ini.AddEntry("VOLUME_CENTER_X",volumeCenter.x);
    ini.AddEntry("VOLUME_CENTER_Y",volumeCenter.y);
    ini.AddEntry("VOLUME_PITCH_X",volumePitch.x);
    ini.AddEntry("VOLUME_PITCH_Y",volumePitch.y);
    ini.ParseFile(filepath);

    std::string btype_str = std::string(btype);
    if( btype_str == std::string("PARALLEL_BEAM") ){
      beamType =BeamType::PARALLEL_BEAM;
    }else if( btype_str == std::string("CONE_BEAM") ){
      beamType =BeamType::CONE_BEAM;
    }

    std::string dtype_str = std::string(dtype);
    if( dtype_str == std::string("HORIZONTAL") ){
      detectorType =DetectorType::HORIZONTAL;
    }

    sourceOrigin.x =gantryCenter.x;
    sourceOrigin.y =gantryCenter.y + sod;

    detectorCenterPos.x =gantryCenter.x;
    detectorCenterPos.y =gantryCenter.y + sod - sid;

  }


  void Print(){
    std::cout << "SID:" << sid <<std::endl;
    std::cout << "SOD:" << sod <<std::endl;

    std::cout << "BEAM_TYPE:";
    switch(beamType){
    case BeamType::PARALLEL_BEAM:
      std::cout << "PARALLEL_BEAM" <<std::endl;
      break;
    case BeamType::CONE_BEAM:
      std::cout << "CONE_BEAM" <<std::endl;
      break;
    default:
      std::cout << "UNKNOWN" << std::endl;
      break;
    }

    std::cout << "GANTRY_CENTER_X:" << gantryCenter.x << std::endl;
    std::cout << "GANTRY_CENTER_Y:" << gantryCenter.y << std::endl;

    std::cout << "DETECTOR_TYPE:";
    switch(detectorType){
    case DetectorType::HORIZONTAL:
      std::cout << "HORIZONTAL" <<std::endl;
      break;
    default:
      std::cout << "UNKNOWN" << std::endl;
      break;
    }

    std::cout << "DETECTOR_CENTER_PIXEL:" << detectorCenterPixel << std::endl;
    std::cout << "DETECTOR_PIXEL_NUM:" << detectorPixelNum << std::endl;
    std::cout << "DETECTOR_PITCH:" << detectorPitch << std::endl;
    std::cout << "PROJECTION_COUNT:" << projectionCount << std::endl;

    std::cout << "VOLUME_SIZE_X:" << volumeSize.x << std::endl;
    std::cout << "VOLUME_SIZE_Y:" << volumeSize.y << std::endl;
    std::cout << "VOLUME_CENTER_X:" << volumeCenter.x << std::endl;
    std::cout << "VOLUME_CENTER_Y:" << volumeCenter.y << std::endl;
    std::cout << "VOLUME_PITCH_X:" << volumePitch.x << std::endl;
    std::cout << "VOLUME_PITCH_Y:" << volumePitch.y << std::endl;
  }

  Real2_t getSourceOrigin(int detectorIndex) const {
    if(beamType == BeamType::PARALLEL_BEAM){
      Real2_t ret =this->getDetectorPos(detectorPixelNum);
      ret.y =ret.y +sid;
      return ret;
    }else{
      return sourceOrigin;
    }
  }


  Real2_t getDetectorPos(int detectorIndex) const {
    if(detectorType == DetectorType::HORIZONTAL){
      Real2_t ret =this->detectorCenterPos;
      ret.x = detectorPitch * (detectorIndex - this->detectorPixelNum / 2 + 0.5F);
      ret.y = gantryCenter.y + sod - sid;
      return ret;
    }else{
      assert(!"invalid detector type");
      return detectorCenterPos;
    }
  }

};
