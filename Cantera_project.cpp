#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include <iostream>
#include <cmath>

using namespace Cantera;

double crank_angle_rad(double t, double engine_speed) 
{
    return engine_speed * t;      // Crank angle [rad]
}

double crank_angle_degree(double t, double engine_speed)
{
    return engine_speed * 180 / Pi * t;      // Crank angle [deg]
}

double piston_speed(double (*ptr), double engine_speed, double piston_stroke)
{       
    return piston_stroke * engine_speed * sin(*ptr);           // Approximate piston speed with sinusoidal velocity profile
}
int main()
{
    addDirectory("C:/Users/JA/anaconda3/envs/ct-dev/share/cantera/data");

        // Input Parameters

        // reaction mechanism, kinetics type and compositions
    std::string reaction_mechanism = "nDodecane_Reitz.yaml";
    std::string phase_name = "nDodecane_IG";
    std::string comp_air = "o2:1, n2 : 3.76";
    std::string comp_fuel = "c12h26:1";

        // Engine Parameters

    double engine_speed = 3000 * 2 * Pi / 60;                                // engine speed[rad / s](3000 rpm)
    double piston_volume = 5 * pow(10, -4);                                  // displaced volume[m * *3]
    double epsilon = 21;                                                     // compression ratio[-]
    double d_piston = 0.083;                                                 // piston diameter[m]
    double piston_stroke = piston_volume / (Pi * pow(d_piston, 2) / 4);      //piston stroke [m]

    //simulation:

    double dt = 1 / (360 * engine_speed / (2 * Pi));
    double sim_time = 0;
    while (sim_time < 360 * dt)
    {



        sim_time += dt;
    }


        return 0;
}
