//================================================//
//					  GPS Tools					  //
//================================================//

// This class provides some utilities in dealing with
// GPS data, by converting between coordinate types
// and finding distances bewteen points
//
// 2015 Enrico Eberhard

#include "GPSTools.h"

#define UTM_a       6378137.0f      //equatorial radius
#define UTM_k0      0.9996f         //scale factor
#define UTM_b       6356752.3142f   //polar radius
#define UTM_e       0.08182f        //ellipse UTM_eentricity sqrt(1-UTM_b^2/UTM_a^2)
#define UTM_e1sq    0.00674f        //UTM_e^2/(1-UTM_e^2)
#define UTM_n       0.00168f        //(UTM_a-UTM_b)/(UTM_a+UTM_b)

#define USGS_A 0.99832426530f
#define USGS_B 0.00251465688f
#define USGS_C 0.00000263915f
#define USGS_D 0.00000000341f

float GPSTools::distAvg(float lat, float lon)
{
    int idx = 0;
    //shift buffer with new points
    for(idx = 1; idx < GPS_BUFF_L; idx++) {
        gpsBuff[idx - 1] = gpsBuff[idx];
    }
    
    gpsBuff[GPS_BUFF_L - 1].lat = lat;
    gpsBuff[GPS_BUFF_L - 1].lon = lon;
    
    
    //find distance between points
    float distBuff = 0.0f;
    for(idx = 1; idx < GPS_BUFF_L; idx++) {
        distBuff += distBetweenLatLon(gpsBuff[idx - 1], gpsBuff[idx]);
    }
    
    return distBuff / (GPS_BUFF_L - 1);
    
}

void GPSTools::setBuff(float lat, float lon)
{
    for(int idx=0; idx < GPS_BUFF_L; idx++) {
        gpsBuff[idx].lat = lat;
        gpsBuff[idx].lon = lon;
    }
}


void GPSTools::setOrigin(float lat, float lon)
{
    origin.lat = lat;
    origin.lon = lon;
}

float GPSTools::getOrigin(int lat1lon2)
{
    float retval = 404.0f;
    if(lat1lon2 == 1) {
        retval = origin.lat;
    }
    else if (lat1lon2 == 2) {
        retval = origin.lon;
    }
    
    return retval;
}

float GPSTools::distBetweenLatLon(gpsPoint p1, gpsPoint p2)
{
    //Simple haversine function to approximate distance m between
    //two lat/long points as arc over sphere. Best for small distances
    //cA = 2asin(sqrt(sin^2(dLat/2)+cos(lat1)cos(lat2)sin^2(dLong/2)))
    //						(a		+		c		   *       b)
    //d = r*cA
    
    double dLat = ((double)p2.lat - (double)p1.lat) * DEG_TO_RAD; //radians between latitudes
    double dLon = ((double)p2.lon - (double)p1.lon) * DEG_TO_RAD; //radians between longitudes
    
    //equation components
    double a = sin(dLat / 2.0f); a = a * a;
    double b = sin(dLon / 2.0f); b = b * b;
    double c = cos(p1.lat * DEG_TO_RAD) * cos(p2.lat * DEG_TO_RAD);
    
    double centralAngle = 2.0f * asin(sqrt(a + (c * b)));
    
    return (float)(UTM_a * centralAngle); //distance in meters
}
float GPSTools::distBetweenLatLon(float lat1, float lon1, float lat2, float lon2)
{
    gpsPoint p1(lat1, lon1);
    gpsPoint p2(lat2, lon2);
    
    return GPSTools::distBetweenLatLon(p1, p2);
}

float GPSTools::lon2x(float lon)
{
    double dLon = ((double)lon - (double)origin.lon) * DEG_TO_RAD; //radians between longitudes
    
    return (float)(UTM_a * dLon);
}
float GPSTools::lat2y(float lat)
{
    double dLat = ((double)lat - (double)origin.lat) * DEG_TO_RAD; //radians between longitudes
    
    return (float)(UTM_a * dLat);
}

void GPSTools::latlontoXY(float lat, float lon, float *x, float *y)
{
    *x = lon2x(lon);
    *y = lat2y(lat);
}
float GPSTools::x2lon(float x)
{
    float cA = x / (float)UTM_a;
    
    return origin.lon + (cA / DEG_TO_RAD);
}
float GPSTools::y2lat(float y)
{
    float cA = y / (float)UTM_a;
    
    return origin.lat + (cA / DEG_TO_RAD);
}


float GPSTools::angleBetweenLatLon(gpsPoint p1, gpsPoint p2)
{
    //very naive direction vector between points
    //N = 0, E = pi/2, etc
    
    double dLat = ((double)p2.lat - (double)p1.lat) * DEG_TO_RAD; //radians between latitudes
    double dLon = ((double)p2.lon - (double)p1.lon) * DEG_TO_RAD; //radians between longitudes
    
    float x = (float)(UTM_a * dLat);  //x axis points north
    float y = (float)(-UTM_a * dLon); //y axis points east
    return atan2(y,x);
}
float GPSTools::angleBetweenLatLon(float lat1, float lon1, float lat2, float lon2)
{
    gpsPoint p1(lat1, lon1);
    gpsPoint p2(lat2, lon2);
    
    return GPSTools::angleBetweenLatLon(p1, p2);
}


void GPSTools::LatLonToUTM(float latCrd, float lonCrd, float *x, float *y)
{
    //convert a set of coordinates to UTM standard x and y.
    //use for better surface distance/direction
    //using WGS84 datum
    
    double lat = (double)latCrd * DEG_TO_RAD;
    double lon = (double)lonCrd * DEG_TO_RAD;
    
    double mem = 1 - pow(UTM_e * sin(lat), 2);
    
    //radius of curvature in the meridian plane (not used)
    //double rho = UTM_a*(1-UTM_e*UTM_e)/pow(mem,1.5);
    
    double nu = UTM_a / pow(mem, 0.5f); //rad. of curv. perp. to mer. plane
    
    //find zone from lon
    int zone = (int)((180.0f + (lonCrd / 100.0f))/6.0f) + 1;
    
    double long0;
    if ((zone < 31) || (zone > 37)) {
        long0 = (double)((zone*6) - 183)*DEG_TO_RAD;
    }
    else {
        //TODO: special cases for norway and svalbard
        //check latitude to see if in non-standard zones
        //set central meridian accordingly
        long0 = (double)((zone*6) - 183)*DEG_TO_RAD; //proceed for now
    }
    
    double p = lon - long0;
    
    double M = UTM_a * ((USGS_A * lat) - (USGS_B * sin(2.0f * lat))
                      + (USGS_C * sin(4.0f * lat)) - (USGS_D * sin(6.0f * lat)));
    
    mem = UTM_k0 * nu;
    double K1 = M * UTM_k0;
    double K2 = mem * sin(2.0f * lat) / 4.0f;
    double K3 = (mem * sin(lat) * pow(cos(lat), 3) / 24.0f)
        * ((5.0f - pow(tan(lat), 2) + 9 * UTM_e1sq * pow(cos(lat), 2)
        + 4.0f * UTM_e1sq * UTM_e1sq * pow(cos(lat), 4)));
    
    double K4 = mem * cos(lat);
    double K5 = (mem * pow(cos(lat), 3) / 6.0f)
    * (1 - pow(tan(lat), 2) + UTM_e1sq * pow(cos(lat), 2));
    
    *y = (float)(K1 + (K2 * p * p) + (K3 * pow(p,4)));
    *x = (float)((K4 * p) + (K5 * pow(p,3)) + 500000.0f);
}


GPSTools::GPSTools(void)
{
    
}

GPSTools::GPSTools(float lat, float lon)
{
    setOrigin(lat,lon);
}
