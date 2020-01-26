//================================================//
//      Pose Estimator through Inertial Update    //
//            and Monte Carlo resampling          //
//================================================//

// This class estimates the position of the sailboat
// using a Monte Carlo particle filter.
// Particles contain state variables and are updated
// according to a system models.
// The resampling phase stochastically selects particles
// based on sensor readings.
//
// 2015 Enrico Eberhard

#ifndef POSE_EST_H
#define POSE_EST_H

#define N_POINTS    500.0f
#define BOAT_MASS   15.0f
#define Half_R_A    0.30625.0f //0.5*wingArea(0.5 m3)*airDensity(1.225 k/m3)
#define RUD_CNTR    1550.0f

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "GPSTools.h"

#include <random>

class PoseEstimator
{
private:
    
    
    struct pointState				//the data within a single particle
    {
        float x, y, theta, s, b;
    };
    pointState state[N_POINTS];		//set of current particles
    pointState newState[N_POINTS];	//set of future particles
    
    pointState estimate;			//estimate of current state
    
    float getLift(int AoA, float VAir);	//calculate lift force
    float getDrag(int AoA, float VAir);	//calculate drag force
    
    float localWind(float gW, float sA, int i); //wind direction relative to boat
    
    //randomly distribute points in space
    float scatter(void);
    
    //with defined variance
    float scatter(float s);
    
    //with mean and variance for x and y coordinates
    void scatter(float *x, float sx, float *y, float sy);
    
    
    
public:
    
    //get global wind direction
    float globalWind(float wA, float sA);
    
    //no argument, get estimate of:
    float getX(void);		//X coord
    float getY(void);		//Y coord
    float getAngle(void);	//heading
    float getSpeed(void);	//speed
    float getBelief(void);	//certainty of estimate
    
    //with argument, get from nth point:
    float getX(int n);		//X coord
    float getY(int n);		//Y coord
    float getAngle(int n);	//heading
    float getSpeed(int n);	//speed
    float getBelief(int n);	//certainty of point
    
    //set up particles with initial heuristic
    void initializeParticles(float x, float y, float yaw);
    
    //update the state of particles from stochastically modelled dynamics
    void updateState(int rudA, float dt);
    
    //resample particles probabilistically based on belief and sensor data
    void resampleParticles(float x, float y, float HDoP, float IMU);
    
    //constructor
    PoseEstimator(void);
    //constructor + initialize particles
    PoseEstimator(float x, float y, float yaw);
    
    
};

#endif
