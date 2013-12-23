#pragma once

// This is redundant with the above but whatever
#ifndef _FLUIDSIM2D_H
#define _FLUIDSIM2D_H

#include <SFML\Graphics.hpp>

// Component enum for better readbility re: boundary constraints
enum Component {
    NONE,
    VERT,
    HORIZ,
    DEPTH
};

// 2D Fluid sim implementation
class FluidSim2D
{
public:
    FluidSim2D(int width, int height, float dt, float diff, float visc);
    ~FluidSim2D();

// Private vars
private:
    int width, height;
    float dt, diff, visc;

    // Velocity arrays
    float *u_l, *v_l, *uPrev, *vPrev;

    // Density arrays
    float *dens_l, *densPrev;

// Public functions
public:
    // A single simulation step
    void simulate();

    // Draws the fluid simulation results out to the window
    void draw(int winX, int winY);

    // Handles UI Input, currently triggered on left click
    void setSource(sf::Vector2i mousePos, int winX, int winY, float source = 100.0f);

    // Sets velocities
    void setVelocity(sf::Vector2i mousePos, int winX, int winY, float force = 5.0f);

    // Resets the source arrays
    void reinit();

// Private functions
private:
    // Density simulation step
    void densStep(float *x, float *x0, float *u, float *v);

    // Velocity simulation step
    void velStep(float *u, float *v, float *u0, float *v0);

    // Adds the source in s to the values in x
    void addSource(float *x, float *s);

    // Diffuses the fluid
    void diffuse(Component b, float *x, float *x0, float d);

    // Advects the fluid
    void advect(Component b, float *d, float *d0, float *u, float *v);

    // Mass-conservation projection function thing
    // AKA: the black magic function
    void project(float *u, float *v, float *p, float *div);

    // Handles boundary restrictions
    void setBound(Component b, float *x);

    void drawDensity();

    void drawVelocity();

    // This is where I wish C++ had properties like C++ does
    int size();
    float hSpacing();
    float vSpacing();

    void debugPrintGrid(float *g);
};

#endif