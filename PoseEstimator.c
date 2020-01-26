//================================================//
//      Pose Estimator through Inertial Update    //
//            and Monte Carlo resampling          //
//================================================//

// This class estimates the position of the the boat
// using a Monte Carlo particle filter.
// Particles contain state variables and are updated
// according to a system models.
// The resampling phase stochastically selects particles
// based on sensor readings.
//
// 2015 Enrico Eberhard

#include "PoseEstimator.h"

PoseEstimator::PoseEstimator(void)
{
    timespec timeStamp;
    clock_gettime(CLOCK_MONOTONIC, &timeStamp);
    srand((uint)timeStamp.tv_nsec);
}

PoseEstimator::PoseEstimator(float x, float y, float yaw)
{
    initializeParticles(x, y, yaw);
}



void PoseEstimator::initializeParticles(float x, float y, float yaw)
{
    float startProb = 10000.0f / (float)N_POINTS;
    
    float sumX = 0.0f;
    float sumY = 0.0f;
    float sumThetaX = 0.0f;
    float sumThetaY = 0.0f;
    float sumS = 0.0f;
    
    
    for (int idx = 0; idx < N_POINTS; idx++) {
        scatter(&state[idx].x, 20, &state[idx].y, 20);
        state[idx].x += x;
        state[idx].y += y;
        state[idx].theta = scatter(30) + yaw;
        state[idx].s = scatter(1.5) + 1;
        state[idx].b = startProb;
        
        while(state[idx].theta >= 180.0f) {
            state[idx].theta -= 360.0f;
        }
        while(state[idx].theta < -180.0f) {
            state[idx].theta += 360.0f;
        }
        
        sumX += state[idx].x;
        sumY += state[idx].y;
        sumThetaX += cos(state[idx].theta * DEG_TO_RAD);
        sumThetaY += sin(state[idx].theta * DEG_TO_RAD);
        sumS += state[idx].s;
    }
    
    estimate.x = sumX / (float)N_POINTS;
    estimate.y = sumY / (float)N_POINTS;
    estimate.theta = atan2((sumThetaY / (float)N_POINTS), (sumThetaX / (float)N_POINTS)) / DEG_TO_RAD;
    estimate.s = sumS / (float)N_POINTS;
}


void PoseEstimator::updateState(int rudA, float dt)
{
    float turnDT = ((float)rudA - RUD_CNTR) / 40.0f;
    
    float wsumX = 0.0f;
    float wsumY = 0.0f;
    float wsumThetaX = 0.0f;
    float wsumThetaY = 0.0f;
    float wsumS = 0.0f;
    float wsumB = 0.0f;
    
    for(int count = 0; count < N_POINTS; count++) {
        state[count].theta = state[count].theta - (turnDT * state[count].s);
        
        state[count].x = state[count].x + (state[count].s * dt * sin(state[count].theta * DEG_T0_RAD));
        state[count].y = state[count].y + (state[count].s * dt * cos(state[count].theta * DEG_T0_RAD));
        
        wsumX += state[count].x * state[count].b;
        wsumY += state[count].y * state[count].b;
        wsumThetaX += cos(state[count].theta * DEG_T0_RAD) * state[count].b;
        wsumThetaY += sin(state[count].theta * DEG_T0_RAD) * state[count].b;
        wsumS += state[count].s * state[count].b;
        wsumB += state[count].b;
    }
    
    //get weighted average
    estimate.x = wsumX / wsumB;
    estimate.y = wsumY / wsumB;
    estimate.theta = atan2(wsumThetaY / wsumB, wsumThetaX / wsumB) / DEG_TO_RAD;
    estimate.s = wsumS / wsumB;
}

void PoseEstimator::resampleParticles(float x, float y, float HDoP, float IMU)
{
    static int pos_scale = 1;
    static float rot_scale = 1.0f;
    
    static float BeliefRatio = 0.8f;
    
    //weighted
    float wsumX = 0.0f;
    float wsumY = 0.0f;
    float wsumThetaX = 0.0f;
    float wsumThetaY = 0.0f;
    float wsumS = 0.0f;
    float wsumB = 0.0f;
    
    //update beliefs
    
    float agrBelief[N_POINTS];
    float sumBelief = 0.0f;
    float newSumBelief = 0.0f;
    
    for(int counter = 0; counter < N_POINTS; counter++) {
        
        //------ GPS --------
        float xDiff = state[counter].x - x;
        float yDiff = state[counter].y - y;
        float offset = sqrt((xDiff * xDiff) + (yDiff * yDiff));
        
        float PosBelief = (float)pos_scale * HDoP / offset;
        
        
        //------ IMU --------
        offset = fabs(IMU - state[counter].theta);
        
        while (offset > 180.0f) {
            offset = fabs(offset - 360.0f);
        }
        
        float ThetaBelief = rot_scale / offset;
        
        //update belief
        float prodBelief = PosBelief * PosBelief * ThetaBelief * 10000.0f;
        float newBelief = (prodBelief * BeliefRatio) * ((1 - BeliefRatio) * state[counter].b);
        
        state[counter].b = newBelief; //New belief of particle
        
        sumBelief += newBelief;
        agrBelief[counter]= sumBelief;
    }
    
    //-------Resample-----------------------
    
    float resRatio = round(float(N_POINTS)*0.95f);
    
    int point = 0;
    
    for (point = 0; point < N_POINTS; point++) {
        
        float r = (float)(rand() % (int)((sumBelief + 1.0f) * 10.0f)) / 10.0f;
        
        int idx = 0;
        while (r > agrBelief[idx])
        {
            idx++;
        }
        
        if (point < resRatio)
        {
            //copy indexed point with normal scatter
            newState[point].x = scatter(0.5f) + state[idx].x;
            newState[point].y = scatter(0.5f) + state[idx].y;
            newState[point].theta = scatter(10.0f) + state[idx].theta;
            
            while(newState[point].theta >= 180.0f) {
                newState[point].theta -= 360.0f;
            }
            while(newState[point].theta < -180.0f) {
                newState[point].theta += 360.0f;
            }
            
            newState[point].s = scatter(0.2f) + state[idx].s;
            
            newState[point].b = state[idx].b;
        }
        else
        {
            //scatter rest around GPS
            newState[point].x = scatter(5.0f) + x;
            newState[point].y = scatter(5.0f) + y;
            newState[point].theta = scatter(20.0f) + IMU;
            
            while(newState[point].theta >= 180.0f) {
                newState[point].theta -= 360.0f;
            }
            while(newState[point].theta < -180.0f) {
                newState[point].theta += 360.0f;
            }
            
            newState[point].s = scatter(1.0f);
            newState[point].b = state[idx].b;
        }
        newSumBelief += newState[point].b;
    }
    
    for(point = 0; point < N_POINTS; point++)
    {
        newState[point].b = newState[point].b / newSumBelief;
        
        wsumX += newState[point].x * newState[point].b;
        wsumY += newState[point].y*newState[point].b;
        wsumThetaX += cos(newState[point].theta * DEG_T0_RAD) * newState[point].b;
        wsumThetaY += sin(newState[point].theta * DEG_T0_RAD) * newState[point].b;
        wsumS += newState[point].s * newState[point].b;
        wsumB += newState[point].b;
    }
    
    //get weighted average
    estimate.x = wsumX / wsumB;
    estimate.y = wsumY / wsumB;
    estimate.theta = atan2(wsumThetaY / wsumB, wsumThetaX / wsumB) / DEG_TO_RAD;
    estimate.s = wsumS / wsumB;
    
    memcpy(state, newState, N_POINTS*sizeof(pointState));
}


float PoseEstimator::scatter(float s)
{
    
    float u = (float)((rand() % 1000.0f) + 1.0f) / 1000.0f;
    float v = (float)((rand() % 1000.0f) + 1.0f) / 1000.0f;
    
    float x = sqrt(-2 * log(u)) * cos(M_PI_2 * v);
    
    float r = rand() % 2;
    
    if (r < 1) {
        x = -x;
    }
    
    return x*s;
}

float PoseEstimator::scatter(void)
{
    return scatter(1.0f);
}


void PoseEstimator::scatter(float *x, float sx, float *y, float sy)
{
    float u = (float)((rand() % 1000.0f) + 1.0f) / 1000.0f;
    float v = (float)((rand() % 1000.0f) + 1.0f) / 1000.0f;
    
    *x = sqrt(-2 * log(u)) * cos(M_PI_2 * v) * sx;
    *y = sqrt(-2 * log(u)) * sin(M_PI_2 * v) * sy;
    
    float r = rand() % 2;
    
    if (r < 1) {
        *x = -*x;
    }
    
    r = rand() % 2;
    if (r < 1) {
        *y = -*y;
    }
}


float PoseEstimator::getLift(int AoA, float VAir)
{
    float cL;
    int pos = -1;
    if (AoA > 180.0f) {
        pos = 1;
    }
    
    AoA = abs(AoA - 180.0f);
    
    if (AoA <= 15.0f) {	//linearly increasing model of cL
        cL = AoA / 15.0f;
    }
    else {            //fast linear stall
        cL = abs((20.0f - AoA) / 5.0f);
    }
    
    return Half_R_A * VAir * VAir * cL * pos;
}

float PoseEstimator::getDrag(int AoA, float VAir)
{
    float cD;
    
    AoA = abs(AoA - 180.0f);
    if (AoA <= 90.0f) {	//parabolic model of cL
        cD = 0.01f + 0.75f * pow(((float)AoA / 40.0f), 3);
    }
    else {			//mirror for backwards wing
        cD = 0.01f + 0.75f * pow(((180.0f - (float)AoA) / 40.0f), 3);
        
        if(cD < 0.0f) {
            cD = 0.0f;
        }
    }
    
    return Half_R_A * VAir * VAir * cD;
}

float PoseEstimator::globalWind(float wA, float sA)
{
    float gW = wA + sA - estimate.theta;
    
    while (gW >= 360.0f) {
        gW -= 360.0f;
    }
    while (gW < 0.0f) {
        gW += 360.0f;
    }
    
    return gW;
}

float PoseEstimator::localWind(float gW, float sA, int i)
{
    float wA = state[i].theta + gW - sA;
    
    while (wA >= 360.0f) {
        wA -= 360.0f;
    }
    while (wA < 0.0f) {
        wA += 360.0f;
    }
    
    return wA;
}


float PoseEstimator::getX(void)
{
    return estimate.x;
}

float PoseEstimator::getY(void)
{
    return estimate.y;
}

float PoseEstimator::getAngle(void)
{
    return estimate.theta;
}

float PoseEstimator::getSpeed(void)
{
    return estimate.s;
}

float PoseEstimator::getBelief(void)
{
    return estimate.b;
}


float PoseEstimator::getX(int n)
{
    float retval = 0.0f;
    if ((n >= 0) && (n < N_POINTS)) {
        retval = state[n].x;
    }
    return retval;
}

float PoseEstimator::getY(int n)
{
    float retval = 0.0f;
    if ((n >= 0) && (n < N_POINTS)) {
        retval = state[n].y;
    }
    return retval;
}

float PoseEstimator::getAngle(int n)
{
    float retval = 0.0f;
    if ((n >= 0) && (n < N_POINTS)) {
        retval = state[n].theta;
    }
    return retval;
}

float PoseEstimator::getSpeed(int n)
{
    float retval = 0.0f;
    if ((n >= 0) && (n < N_POINTS)) {
        retval = state[n].s;
    }
    return retval;
}

float PoseEstimator::getBelief(int n)
{
    float retval = 0.0f;
    if ((n >= 0) && (n < N_POINTS)) {
        retval = state[n].b;
    }
    return retval;
}
