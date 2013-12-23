#pragma once

// This is redundant with the above but whatever
#ifndef _FLUIDSIM3D_H
#define _FLUIDSIM3D_H

#include <SFML\Graphics.hpp>
#include <SFML\OpenGL.hpp>
#include <glext.h>
#include "FluidSim2D.h"

// 2D Fluid sim implementation
class FluidSim3D
{
public:
    FluidSim3D(int width, int height, int depth, float dt, float diff, float visc, bool bounds = false);
    ~FluidSim3D();

// Private vars
private:
    int width, height, depth;
    float dt, diff, visc;

    // Velocity arrays
	// _Old are used for vortices, they store the last velocity to calculate the 
	//   gradient, may have been able to use _Prev but added complexity
    float *u_l, *v_l, *w_l, *uPrev, *vPrev, *wPrev, *uOld, *vOld, *wOld;

    // Density arrays
    float *dens_l, *densPrev;

// Public functions
public:
    // Doesn't really do nice object oriented things with these but eh.
    int xrot, yrot, zrot;
    float dist;
    bool bounds; // Turns bounds on or off. Defaults to off.

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
	// adds vortices uses omega...
	void addVortices(float *x, float *y, float *z, float *xPrev, float *yPrev, float *zPrev);

    // Density simulation step
    void densStep(float *x, float *x0, float *u, float *v, float *w);

    // Velocity simulation step
    void velStep(float *u, float *v, float *w, float *u0, float *v0, float *w0);

    // Adds the source in s to the values in x
    void addSource(float *x, float *s);

    // Diffuses the fluid
    void diffuse(Component b, float *x, float *x0, float d);

    // Advects the fluid
    void advect(Component b, float *d, float *d0, float *u, float *v, float *w);

    // Mass-conservation projection function thing
    // AKA: the black magic function
    void project(float *u, float *v, float *w, float *p, float *div);

    // Handles boundary restrictions
    void setBound(Component b, float *x);

    void drawBoundingCube();
    void drawDensity();
    void drawVelocity(int winX, int winY);

    // This is where I wish C++ had properties like C++ does
    int size();
    float hSpacing();
    float vSpacing();
	float dSpacing();

    void debugPrintGrid(float *g);

    // Here comes the GL...
    GLuint fluidTexID;   // Handle for the GL texture ID of the density
    char* RGBABuffer;    // Space for raw data for the densities
    PFNGLTEXIMAGE3DPROC glTexImage3D;   // Black magic function pointer
};

#endif