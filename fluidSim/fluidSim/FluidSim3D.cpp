#include "FluidSim3D.h"
#include <SFML\OpenGL.hpp>
#include <iostream>
#include <glext.h>
#include <math.h>

#define IDX(x,y,z) ((x) + ((width + 2) * (y)) + ((width+2) * (height+2) * (z)))

FluidSim3D::FluidSim3D(int width, int height, int depth, float dt, float diff, float visc, bool bounds) :
    width(width), height(height), depth(depth), dt(dt), diff(diff), visc(visc),
    xrot(0), yrot(0), zrot(0), dist(2.5), bounds(bounds)
{
    u_l = new float[size()]();
    v_l = new float[size()]();
    w_l = new float[size()]();

    uPrev = new float[size()]();
    vPrev = new float[size()]();
    wPrev = new float[size()]();

    uOld = new float[size()]();
    vOld = new float[size()]();
    wOld = new float[size()]();

    dens_l = new float[size()]();
    densPrev = new float[size()]();

    // Generate texture ID for visulalization
    glGenTextures(1, &fluidTexID);

    // Allocate a buffer for it to get the updated densities
    RGBABuffer = new char[size() * 4];

    // Apparently glTexImage3D isn't standard so we do some black magic and load it from
    // a magic GL .dll somewhere. It didn't crash when I ran it so I guess it works?
    glTexImage3D = (PFNGLTEXIMAGE3DPROC) wglGetProcAddress("glTexImage3D");
}


FluidSim3D::~FluidSim3D(void)
{
    delete u_l;
    delete v_l;
	delete w_l;
    delete uPrev;
    delete vPrev;
	delete wPrev;
    delete dens_l;
    delete densPrev;
}

// Runs a simulation step
void FluidSim3D::simulate()
{
    // Updates from UI come in in the main SFM loop
	velStep(u_l, v_l, w_l, uPrev, vPrev, wPrev);
	
	addVortices(u_l, v_l, w_l, uOld, vOld, wOld);
    
	densStep(dens_l, densPrev, u_l, v_l, w_l);
}

void FluidSim3D::draw(int winX, int winY)
{
    // Happy fun time drawing code here
    glViewport(0, 0, winX, winY);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // I can't recommend velocities and densities being drawn at the same time
    // Densities use an ortho projection and velocities use a perspective
    drawDensity();

	//drawVelocity(winX, winY);

    glFlush();
}

void FluidSim3D::drawBoundingCube()
{
    glBegin(GL_LINES);
        glColor3f(0.0, 0.0, 1.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 1.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 1.0, 1.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, 1.0);
        glVertex3f(0.0, 0.0, 1.0);
        glVertex3f(1.0, 0.0, 1.0);
        glVertex3f(0.0, 0.0, 1.0);
        glVertex3f(0.0, 1.0, 1.0);
        glVertex3f(1.0, 1.0, 1.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(1.0, 1.0, 1.0);
        glVertex3f(0.0, 1.0, 1.0);
        glVertex3f(1.0, 1.0, 1.0);
        glVertex3f(1.0, 0.0, 1.0);
    glEnd();
}

void FluidSim3D::drawDensity()
{
    glOrtho(-1.0, 1.0, -1.0, 1.0, -2.0, 2.0);
    glMatrixMode(GL_PROJECTION);

    // Derived from http://www.codeproject.com/Articles/352270/Getting-started-with-Volume-Rendering
    // And from http://physbam.stanford.edu/~fedkiw/papers/stanford2001-01.pdf

    // Convert densities to something OpenGL understands
    glBindTexture(GL_TEXTURE_3D, fluidTexID);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    // For simplicity, we're going to assume that the light source in the scene is directly above the volume
    // Compute transmission based on a light source to simulate self-shadowing

    // So I don't actually know how the extinction cross section coefficient is calculated, but I'll just pick a value here...
    double extc = 1.0;
    double albedo = 0.98;
    double lLight = 10.0;

    for (int i = 0; i < width + 2; i++)
    {
        for (int k = 0; k < depth + 2; k++)
        {
            double tRay = 1.0;

            for (int j = height + 1; j >= 0; j--)
            {
                double cext;
                
                if (dens_l[IDX(i, j, k)] == 0)
                    cext = 0;
                else
                    cext = extc * dens_l[IDX(i, j, k)];

                double tVox = exp(-cext * vSpacing());
                double lVox = (albedo * lLight * (1.0 - tVox)) * tRay;

                tRay *= tVox;

                lVox = (lVox > 1.0) ? 1.0 : lVox;
                tVox = (tVox > 1.0) ? 1.0 : tVox;

                char color = (char) (lVox * 255);
                char opacity = (char) ((1 - tVox) * 255);
                RGBABuffer[IDX(i, j, k) * 4] = color;
                RGBABuffer[IDX(i, j, k) * 4 + 1] = color;
                RGBABuffer[IDX(i, j, k) * 4 + 2] = color;
                RGBABuffer[IDX(i, j, k) * 4 + 3] =  opacity;
            }
        }
    }

    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, width + 2, height + 2, depth + 2, 0, GL_RGBA, GL_UNSIGNED_BYTE, RGBABuffer);
    glBindTexture(GL_TEXTURE_3D, 0);

    // Render party time
    glEnable(GL_ALPHA_TEST);
    glAlphaFunc(GL_GREATER, 0.03f);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glMatrixMode(GL_TEXTURE);
    glLoadIdentity();

    glTranslatef(0.5, 0.5, 0.5);
    glRotatef(xrot, 1, 0, 0);
    glRotatef(yrot, 0, 1, 0);
    glRotatef(zrot, 0, 0, 1);
    glTranslatef(-0.5, -0.5, -0.5);

    glEnable(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, fluidTexID);

    glBegin(GL_QUADS);
        for (float ind = 1.0f; ind >= -1.0f; ind -= 0.005f)
        {
            glTexCoord3f(-0.15, -0.15, (ind + 1.0f) / 2.0f);
            glVertex3f(-1, -1, ind);
            glTexCoord3f(1.15, -0.15, (ind + 1.0f) / 2.0f);
            glVertex3f(1, -1, ind);
            glTexCoord3f(1.15, 1.15, (ind + 1.0f) / 2.0f);
            glVertex3f(1, 1, ind);
            glTexCoord3f(-0.15, 1.15, (ind + 1.0f) / 2.0f);
            glVertex3f(-1, 1, ind);
        }
    glEnd();
}

void FluidSim3D::drawVelocity(int winX, int winY)
{
    gluPerspective(60, (GLdouble) winX / (GLdouble) winY, 0.1, 1000);
    gluLookAt(0.5, 0.5, dist, 0.5, 0.5, 0.5, 0, 1, 0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0.5, 0.5, 0.5);
    glRotatef(xrot, 1, 0, 0);
    glRotatef(yrot, 0, 1, 0);
    glRotatef(zrot, 0, 0, 1);
    glTranslatef(-0.5, -0.5, -0.5);

    drawBoundingCube();

	int i, j, k;
	float x, y, z, hx, hy, hz;

	hx = hSpacing();
    hy = vSpacing();
	hz = dSpacing();

	glColor3f(1.0f, 0, 0);
	glLineWidth (1.0f);

	glBegin (GL_LINES);
		for (i = 1; i <= width ; i++ )
        {
			x = (i - 0.5f) * hx;
			for (j = 1; j <= height; j++ )
            {
				y = (j - 0.5f) * hy;
				for (k = 1; k <= height; k++ )
				{
					z = (k - 0.5f) * hz;

					glVertex3f(x,y,z);
					glVertex3f(x + u_l[IDX(i,j,k)], y + v_l[IDX(i,j,k)], z + w_l[IDX(i,j,k)]);
				}
			}
		}
	glEnd();
}

void FluidSim3D::setSource(sf::Vector2i mousePos, int winX, int winY, float source)
{
    int x = (int)((mousePos.x / (float)winX) * width);
    int y = (int)(((winY - mousePos.y) / (float)winY) * height);

    if (x < 1 || x > width || y < 1 || y > height)
        return;

	// This will put all sources on the first plane
    densPrev[IDX(x,y,2)] = source;
}

void FluidSim3D::setVelocity(sf::Vector2i mousePos, int winX, int winY, float force)
{
    int x = (int)((mousePos.x / (float)winX) * width);
    int y = (int)(((winY - mousePos.y) / (float)winY) * height);

    if (x < 1 || x > width || y < 1 || y > height)
        return;

	//upward force/velocity
    wPrev[IDX(x,y,2)] = force;
}

void FluidSim3D::reinit()
{
    for (int i = 0; i < size(); i++)
    {
        densPrev[i] = uPrev[i] = vPrev[i] = wPrev[i] = 0.0f;
    }
}

// will over write values in _Prev with values of omega (save space)
// calling the _Prev is leftover from a previous attempt at this, no actual data is useful in them at begginning or end
void FluidSim3D::addVortices(float *x, float *y, float *z, float *xPrev, float *yPrev, float *zPrev)
{
    int i, j, k;
	float delX, delY, delZ;
	float omegaX, omegaY, omegaZ;
	float NX, NY, NZ;
	float smoke_epsilon = 30;

	float *omegaLen = new float[size()];

	// calculate omega
    for (i = 1; i <= width; i++)
        for (j = 1; j <= height; j++)
			for (k = 1; k <= depth; k++)
			{
				delX = (x[IDX(i+1,j,k)] - x[IDX(i-1,j,k)]) / (hSpacing() * 2);
				delY = (y[IDX(i,j+1,k)] - y[IDX(i,j-1,k)]) / (vSpacing() * 2);
				delZ = (z[IDX(i,j,k+1)] - z[IDX(i,j,k-1)]) / (dSpacing() * 2);
				
				// calculate omega by taking grad(velocity) x velocity
				xPrev[IDX(i,j,k)] = delY*z[IDX(i,j,k)] - delZ*y[IDX(i,j,k)];
				yPrev[IDX(i,j,k)] = delZ*x[IDX(i,j,k)] - delX*z[IDX(i,j,k)];
				zPrev[IDX(i,j,k)] = delX*y[IDX(i,j,k)] - delY*x[IDX(i,j,k)];

				omegaLen[IDX(i,j,k)] = sqrt(xPrev[IDX(i,j,k)]*xPrev[IDX(i,j,k)] 
					+ yPrev[IDX(i,j,k)] * yPrev[IDX(i,j,k)]
					+ zPrev[IDX(i,j,k)] * zPrev[IDX(i,j,k)]);
			}
			
	// now compute N's and forces
    for (i = 1; i <= width; i++)
        for (j = 1; j <= height; j++)
			for (k = 1; k <= depth; k++)
			{
				omegaX = -xPrev[IDX(i,j,k)];
				omegaY = -yPrev[IDX(i,j,k)];
				omegaZ = -zPrev[IDX(i,j,k)];

				// N = normalize(grad(|omega|))
				NX = (omegaLen[IDX(i+1,j,k)] - omegaLen[IDX(i-1,j,k)]) / (hSpacing() * 2);
				NY = (omegaLen[IDX(i,j+1,k)] - omegaLen[IDX(i,j-1,k)]) / (hSpacing() * 2);
				NZ = (omegaLen[IDX(i,j,k+1)] - omegaLen[IDX(i,j,k-1)]) / (hSpacing() * 2);

				float len = sqrt(NX*NX + NY*NY + NZ*NZ);
				if(len==0) {
					NX = 0;
					NY = 0;
					NZ = 0;
				} else {
					NX /= len;
					NY /= len;
					NZ /= len;
				}

				// compute and add force by epsilon * h * (N x omega)
				x[IDX(i,j,k)] += dt * smoke_epsilon * hSpacing() * (NY*omegaZ - NZ*omegaY);
				y[IDX(i,j,k)] += dt * smoke_epsilon * vSpacing() * (NZ*omegaX - NX*omegaZ);
				z[IDX(i,j,k)] += dt * smoke_epsilon * dSpacing() * (NX*omegaY - NY*omegaX);
			}
	
}

// Updates the density.
// Gets sources, diffuses and advects the density
void FluidSim3D::densStep(float *x, float *x0, float *u, float *v, float *w)
{
   addSource(x, x0);
   diffuse(NONE, x0, x, diff);
   advect(NONE, x, x0, u, v, w);
}

// Updates the velocity
void FluidSim3D::velStep(float *u, float *v, float *w, float *u0, float *v0, float *w0)
{
    addSource(u, u0);
    addSource(v, v0);
    addSource(w, w0);

    diffuse(HORIZ, u0, u, visc);
    diffuse(VERT, v0, v, visc);
    diffuse(DEPTH, w0, w, visc);

    project(u0, v0, w0, u, v);

    advect(HORIZ, u, u0, u0, v0, w0);
    advect(VERT, v, v0, u0, v0, w0);
    advect(DEPTH, w, w0, u0, v0, w0);
    
    project(u, v, w, u0, v0);
}

// Adds a source to an existing array of values
void FluidSim3D::addSource(float *x, float *s)
{
    int i, size = this->size();

    for (i = 0; i < size; i++)
        x[i] += dt * s[i];
}

// Solves a linear system to diffuse the fluid. Solution method taken from Stamm paper.
// Possible Change: Different solver?
void FluidSim3D::diffuse(Component b, float *x, float *x0, float d)
{
    int i, j, k, c;
    float a = dt * d * width * height;

    for (c = 0; c < 20; c++)
    {
        for (i = 1; i <= width; i++)
        {
            for (j = 1; j <= height; j++)
            {
				for (k = 1; k <= depth; k++)
				{
					x[IDX(i, j, k)] = (x0[IDX(i, j, k)] + a * (x[IDX(i-1, j, k)] + x[IDX(i+1, j, k)] + x[IDX(i, j-1, k)] + x[IDX(i, j+1, k)] + x[IDX(i, j, k-1)] + x[IDX(i, j, k+1)]))
									/ (1 + 6 * a);
				}
            }
        }
        setBound(b, x);
    }
}

// Does an advection step on the given components
void FluidSim3D::advect(Component b, float *d, float *d0, float *u, float *v, float *w)
{
    int i, j, k, i0, j0, k0, i1, j1, k1;

    float x, y, z, s0, t0, r0, s1, t1, r1, dtx, dty, dtz; // r is for z components

    dtx = dt * width;
    dty = dt * height;
    dtz = dt * depth;

    for (i = 1; i <= width; i++)
    {
        for (j = 1; j <= height; j++)
        {
			for (k = 1; k <= depth; k++)
			{
				// Apparently this is a "simple linear backtrace"
				// x,y,z are locations if you were to go backwards dt*N along the velocity
				x = i - dtx * u[IDX(i, j, k)];
				y = j - dty * v[IDX(i, j, k)];
				z = k - dtz * w[IDX(i, j, k)];
            
				// Keep cells in bounds it looks like
				// _0 and _1 for ijk are the grid squares that bound the location
				if (x < 0.5)
					x = 0.5;
				if (x > width + 0.5)
					x = width + 0.5;
				i0 = (int) x;
				i1 = i0 + 1;

				if (y < 0.5)
					y = 0.5;
				if (y > height + 0.5)
					y = height + 0.5;
				j0 = (int) y;
				j1 = j0 + 1;

				if (z < 0.5)
					z = 0.5;
				if (z > depth + 0.5)
					z = depth + 0.5;
				k0 = (int) z;
				k1 = k0 + 1;

				s1 = x - i0; // _1 = difference from starting cell
				s0 = 1 - s1; // _0 = negative + 1 of _0, corresponds to negative velocity out?
				t1 = y - j0;
				t0 = 1 - t1;
				r1 = z - k0;
				r0 = 1 - r1;

				// add the values in d0 of all 8 cells surrounding the caluculated x,y,z, scale by str
				d[IDX(i, j, k)] = s0 * (t0 * (r0 * d0[IDX(i0, j0, k0)] + r1 * d0[IDX(i0, j0, k1)]) + t1 *(r0 * d0[IDX(i0, j1, k0)] + r1 * d0[IDX(i0, j1, k1)])) + 
					s1 * (t0 * (r0 * d0[IDX(i1, j0, k0)] + r1 * d0[IDX(i1, j0, k1)]) + t1 *(r0 * d0[IDX(i1, j1, k0)] + r1 * d0[IDX(i1, j1, k1)]));
							  // s1 * (t0 * d0[IDX(i1, j0)] + t1 * d0[IDX(i1, j1)]);
			}
        }
    }
    setBound(b, d);
}

// Does some sort of magic projection to get a mass-conserving flow
// It should be noted that p and div are just scratch arrays. They are passed in so that
// the fluid sim doesn't have to waste time allocating new memory every time this function is called.
void FluidSim3D::project(float *u, float *v, float *w, float *p, float *div)
{
	int i, j, k, c;

	float hx, hy, hz;

	hx = hSpacing();
	hy = vSpacing();
	hz = dSpacing();

	for (i = 1; i <= width; i++)
	{
		for (j = 1; j <= height; j++)
		{
			for (k = 1; k <= depth; k++)
			{
				div[IDX(i, j, k)] = -0.5 * (hx * u[IDX(i + 1, j, k)] - hx * u[IDX(i - 1, j, k)] +
											hy * v[IDX(i, j + 1, k)] - hy * v[IDX(i, j - 1, k)] +
											hz * w[IDX(i, j, k + 1)] - hz * w[IDX(i, j, k - 1)]);
				p[IDX(i, j, k)] = 0;
			}
		}
	}
	setBound(NONE, div);
	setBound(NONE, p); // why? should all by 0 - it was in the paper don't mess with it :P

	for (c = 0; c < 20; c++)
	{
		for (i = 1; i <= width; i++)
		{
            for (j = 1; j <= height; j++)
            {
                for (k = 1; k <= depth; k++)
                {
                    p[IDX(i, j, k)] = (div[IDX(i, j, k)] + p[IDX(i - 1, j, k)] + p[IDX(i + 1, j, k)] + p[IDX(i, j - 1, k)] + p[IDX(i, j + 1, k)] + 
                                        p[IDX(i, j, k - 1)] + p[IDX(i, j, k + 1)]) / 6;
                }
            }
        }
		setBound(NONE, p);
	}

	for (i = 1; i <= width; i++)
	{
		for (j = 1; j <= height; j++)
		{
			for (k = 1; k <= depth; k++)
			{
				u[IDX(i, j, k)] -= 0.5 * (p[IDX(i + 1, j, k)] - p[IDX(i - 1, j, k)]) / hx;
				v[IDX(i, j, k)] -= 0.5 * (p[IDX(i, j + 1, k)] - p[IDX(i, j - 1, k)]) / hy;
				w[IDX(i, j, k)] -= 0.5 * (p[IDX(i, j, k + 1)] - p[IDX(i, j, k - 1)]) / hz;
			}
		}
	}
	setBound(HORIZ, u);
	setBound(VERT, v);
	setBound(DEPTH, w); // ?
}

// 3D conversion might not be correct, may want to double check
void FluidSim3D::setBound(Component b, float *x)
{
    // Short circuit the function if we want no bounding conditions
    if (!bounds)
        return;

    int i, j, k; //don't need a 3rd variable, but makes code easier to read

	// So VERT, HORIZ, DEPTH are a bit different in 3D... I'll explain as we go.
    // In this next loop, we consider the yz planes, the ones that the x velocity component
    // needs to be contained by.    
    for (i = 1; i <= height; i++)
    {
        for (j = 1; j <= depth; j++)
        {
            // yz plane at x = 0
            // If we're looking to clamp the HORIZ component, invert the bound velocity
            x[IDX(0, i, j)] = (b == HORIZ) ? -x[IDX(1, i, j)] : x[IDX(1, i, j)];
            
            // yz plane at x = width + 2
            x[IDX(width + 1, i, j)] = (b == HORIZ) ? -x[IDX(width, i, j)] : x[IDX(width, i, j)];
        }
    }

    // On to the next planes, the xz planes
    for (i = 1; i <= width; i++)
    {
        for (j = 1; j <= depth; j++)
        {
            // xz plane at y = 0
            x[IDX(i, 0, j)] = (b == VERT) ? -x[IDX(i, 1, j)] : x[IDX(i, 1, j)];

            // xz plane at y = height + 2
            x[IDX(i, height + 1, j)] = (b == VERT) ? -x[IDX(i, height, j)] : x[IDX(i, height, j)];
        }
    }

    // Final planes, the xy planes
    for (i = 1; i <= width; i++)
    {
        for (j = 1; j <= height; j++)
        {
            // xy plane at z = 0
            x[IDX(i, j, 0)] = (b == DEPTH) ? -x[IDX(i, j, 1)] : x[IDX(i, j, 1)];

            // xy plane at z = depth + 2
            x[IDX(i, j, depth + 1)] = (b == DEPTH) ? -x[IDX(i, j, depth)] : x[IDX(i, j, depth)];
        }
    }

	// calculate edges of cube, average of neighboring cells that are on a face
	for(i = 1; i <= width; i++) {
		x[IDX(i, 0, 0)] = 0.5 * (x[IDX(i, 0, 1)] + x[IDX(i, 1, 0)]); 
		x[IDX(i, 0, depth + 1)] = 0.5 * (x[IDX(i, 0, depth)] + x[IDX(i, 1, depth + 1)]); 
		x[IDX(i, height + 1, 0)] = 0.5 * (x[IDX(i, height + 1, 1)] + x[IDX(i, height, 0)]); 
		x[IDX(i, height + 1, depth + 1)] = 0.5 * (x[IDX(i, height + 1, depth)] + x[IDX(i, height, depth + 1)]); 
	}
	for(j = 1; j <= height; j++) {
		x[IDX(0, j, 0)] = 0.5 * (x[IDX(0, j, 1)] + x[IDX(1, j, 0)]); 
		x[IDX(0, j, depth + 1)] = 0.5 * (x[IDX(0, j, depth)] + x[IDX(1, j, depth + 1)]); 
		x[IDX(width + 1, j, 0)] = 0.5 * (x[IDX(width + 1, j, 1)] + x[IDX(width, j, 0)]); 
		x[IDX(width + 1, j, depth + 1)] = 0.5 * (x[IDX(width + 1, j, depth)] + x[IDX(width, j, depth + 1)]); 
	}
	for(k = 1; k <= depth; k++) {
		x[IDX(0, 0, k)] = 0.5 * (x[IDX(0, 1, k)] + x[IDX(1, 0, k)]); 
		x[IDX(0, height + 1, k)] = 0.5 * (x[IDX(0, height, k)] + x[IDX(1, height + 1, k)]); 
		x[IDX(width + 1, 0, k)] = 0.5 * (x[IDX(width + 1, 1, k)] + x[IDX(width, 0, k)]); 
		x[IDX(width + 1, height + 1, k)] = 0.5 * (x[IDX(width + 1, height, k)] + x[IDX(width, height + 1, k)]);
	}

	// NOTE: I don't know if we actually care about the corners but oh well
	// Corner cases:
	x[IDX(0, 0, 0)] = (x[IDX(1, 0, 0)] + x[IDX(0, 1, 0)] + x[IDX(0, 0, 1)]) / 3;
	x[IDX(0, height + 1, 0)] = (x[IDX(1, height + 1, 0)] + x[IDX(0, height, 0)] + x[IDX(0, height + 1, 1)]) / 3;
	x[IDX(0, 0, depth + 1)] = (x[IDX(1, 0, depth + 1)] + x[IDX(0, 1, depth + 1)] + x[IDX(0, 0, depth)]) / 3;
	x[IDX(0, height + 1, depth + 1)] = (x[IDX(1, height + 1, depth + 1)] + x[IDX(0, height, depth + 1)] + x[IDX(0, height + 1, depth)]) / 3;
	
	x[IDX(width + 1, 0, 0)] = (x[IDX(width, 0, 0)] + x[IDX(width + 1, 1, 0)] + x[IDX(width + 1, 0, 1)]) / 3;
	x[IDX(width + 1, height + 1, 0)] = (x[IDX(width, height + 1, 0)] + x[IDX(width + 1, height, 0)] + x[IDX(width + 1, height + 1, 1)]) / 3;
	x[IDX(width + 1, 0, depth + 1)] = (x[IDX(width, 0, depth + 1)] + x[IDX(width + 1, 1, depth + 1)] + x[IDX(width + 1, 0, depth)]) / 3;
	x[IDX(width + 1, height + 1, depth + 1)] = (x[IDX(width, height + 1, depth + 1)] + x[IDX(width + 1, height, depth + 1)] + x[IDX(width + 1, height + 1, depth)]) / 3;
}

// Calculate the size of the grid. Remember that the outside cells are boundary cells.
inline
int FluidSim3D::size()
{
    return (width + 2) * (height + 2) * (depth + 2);
}

inline
float FluidSim3D::hSpacing()
{
    return 1.0 / (float)width;
}

inline
float FluidSim3D::vSpacing()
{
    return 1.0 / (float)height;
}

inline
float FluidSim3D::dSpacing()
{
    return 1.0 / (float)depth;
}

void FluidSim3D::debugPrintGrid(float *g)
{
    for (int i = 0; i < width + 2; i++) {
        for (int j = 0; j < height + 2; j++) {
			for (int k = 0; k < depth + 2; k++) {
				std::cout << g[IDX(i, j, k)] << " ";
			}
			std::cout << std::endl;
		}
        std::cout << std::endl;
    }
    std::cout << std::endl;
}