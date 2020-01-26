#ifndef GPS_TOOLS_H
#define GPS_TOOLS_H

//================================================//
//					  GPS Tools					  //
//================================================//

// This class provides some utilities in dealing with
// GPS data, by converting between coordinate types
// and finding distances bewteen points
//
// 2015 Enrico Eberhard

#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define GPS_BUFF_L 5
#define DEG_TO_RAD (M_PI / 180.0f)

class GPSTools
{
    
private:
    
    struct gpsPoint	//a point with a latitude and longitude
    {
        float lat, lon;
    };
    
    gpsPoint gpsBuff[GPS_BUFF_L]; 	//a time-series of points
    
    gpsPoint origin;				//a local grid origin point
	
	//distance (meters) between two lat/lon points
    float distBetweenLatLon(gpsPoint p1, gpsPoint p2);
	
	//cardinal angle (radians) between two lat/lon points
    float angleBetweenLatLon(gpsPoint p1, gpsPoint p2);
    
public:
	
	//running average of distances between GPS samples
	//(because pos can be a little jumpy, this provides a
	// speed average knowing GPS sample rate)
    float distAvg(float lat, float lon);
	
	//overwrite running average buffer with defined point.
    void setBuff(float lat, float lon);
	
	//set the local grid origin point
    void setOrigin(float lat, float lon);
	
	//return the local grid origin point
    float getOrigin(int lat1lon2);
	
	//in local grid, x axis is east, y axis is north from origin
    float x2lon(float x); //X position to longitude
    float y2lat(float y); //Y position to latitude
	
	//local X,Y grid coordinates from latitude/longitude
    void latlontoXY(float lat, float lon, float *x, float *y);
    
    float lon2x(float lon); //X position from longitude
	float lat2y(float lat); //Y position from latitude
	
	//distance (meters) between two lat/lon points
    float distBetweenLatLon(float lat1, float lon1, float lat2, float lon2);
	
	//cardinal angle (radians) between two lat/lon points
    float angleBetweenLatLon(float lat1, float lon1, float lat2, float lon2);
	
	//UTM X,Y coordinates from latitude/longitude
    void LatLonToUTM(float latCrd, float lonCrd, float *x, float *y);
	
	//constructors
    GPSTools(void);
    GPSTools(float lat, float lon); //sets local origin point
};

#endif
