#include "FluidSim2D.h"
#include <SFML\OpenGL.hpp>
#include <iostream>

#define IDX(x,y) ((x) + ((width + 2) * (y)))

FluidSim2D::FluidSim2D(int width, int height, float dt, float diff, float visc) :
    width(width), height(height), dt(dt), diff(diff), visc(visc)
{
    u_l = new float[size()]();
    v_l = new float[size()]();
    uPrev = new float[size()]();
    vPrev = new float[size()]();
    dens_l = new float[size()]();
    densPrev = new float[size()]();
}


FluidSim2D::~FluidSim2D(void)
{
    delete u_l;
    delete v_l;
    delete uPrev;
    delete vPrev;
    delete dens_l;
    delete densPrev;
}

// Runs a simulation step
void FluidSim2D::simulate()
{
    // Updates from UI come in in the main SFM loop
    velStep(u_l, v_l, uPrev, vPrev);
    densStep(dens_l, densPrev, u_l, v_l);
}

void FluidSim2D::draw(int winX, int winY)
{
    // Happy fun time drawing code here
    glViewport(0, 0, winX, winY);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    drawDensity();
    drawVelocity();

    glFlush();
}

void FluidSim2D::drawDensity()
{
    // Density drawing wooo. Borrowed from Stam's demo because I don't want to spend much
    // time on figureing out rendering right now
    int i, j;
    float x, y, d00, d01, d10, d11;

    float hx = hSpacing();
    float hy = vSpacing();
    
    glBegin(GL_QUADS);

        for (j = 0; j <= height; j++)
        {
            y = (j - 0.5) * hy;
            for (i = 0; i <= width; i++)
            {
                x = (i - 0.5) * hx;
                d00 = dens_l[IDX(i, j)];
                d01 = dens_l[IDX(i, j+1)];
                d10 = dens_l[IDX(i + 1,j)];
                d11 = dens_l[IDX(i + 1, j + 1)];

                glColor3f ( d00, d00, d00 ); glVertex2f ( x, y );
                glColor3f ( d10, d10, d10 ); glVertex2f ( x+hx, y );
                glColor3f ( d11, d11, d11 ); glVertex2f ( x+hx, y+hy );
                glColor3f ( d01, d01, d01 ); glVertex2f ( x, y+hy );
            }
        }

    glEnd();
}

void FluidSim2D::drawVelocity()
{
	int i, j;
	float x, y, hx, hy;

	hx = hSpacing();
    hy = vSpacing();

	glColor3f(1.0f, 0, 0);
	glLineWidth (1.0f);

	glBegin (GL_LINES);
        for (j = 1; j <= height; j++ )
        {
            y = (j - 0.5f) * hy;
            for (i = 1; i <= width ; i++ )
            {
                x = (i - 0.5f) * hx;

				glVertex2f(x,y);
				glVertex2f(x + u_l[IDX(i,j)], y + v_l[IDX(i,j)]);
			}
		}
	glEnd();
}

void FluidSim2D::setSource(sf::Vector2i mousePos, int winX, int winY, float source)
{
    int x = (int)((mousePos.x / (float)winX) * width);
    int y = (int)(((winY - mousePos.y) / (float)winY) * height);

    if (x < 1 || x > width || y < 1 || y > height)
        return;

    densPrev[IDX(x,y)] = source;
}

void FluidSim2D::setVelocity(sf::Vector2i mousePos, int winX, int winY, float force)
{
    int x = (int)((mousePos.x / (float)winX) * width);
    int y = (int)(((winY - mousePos.y) / (float)winY) * height);

    if (x < 1 || x > width || y < 1 || y > height)
        return;

	//upward force/velocity
    vPrev[IDX(x,y)] = force;
}

void FluidSim2D::reinit()
{
    for (int i = 0; i < size(); i++)
    {
        densPrev[i] = uPrev[i] = vPrev[i] = 0.0f;
    }
}

// Updates the density.
// Gets sources, diffuses and advects the density
void FluidSim2D::densStep(float *x, float *x0, float *u, float *v)
{
   addSource(x, x0);
   diffuse(NONE, x0, x, diff);
   advect(NONE, x, x0, u, v);
}

// Updates the velocity
void FluidSim2D::velStep(float *u, float *v, float *u0, float *v0)
{
    addSource(u, u0);
    addSource(v, v0);

    diffuse(VERT, u0, u, visc);
    diffuse(HORIZ, v0, v, visc);

    project(u0, v0, u, v);

    advect(VERT, u, u0, u0, v0);
    advect(HORIZ, v, v0, u0, v0);
    
    project(u, v, u0, v0);
}

// Adds a source to an existing array of values
void FluidSim2D::addSource(float *x, float *s)
{
    int i, size = this->size();

    for (i = 0; i < size; i++)
        x[i] += dt * s[i];
}

// Solves a linear system to diffuse the fluid. Solution method taken from Stamm paper.
// Possible Change: Different solver?
void FluidSim2D::diffuse(Component b, float *x, float *x0, float d)
{
    int i, j, k;
    float a = dt * d * width * height;

    for (k = 0; k < 20; k++)
    {
        for (j = 1; j <= height; j++)
        {
            for (i = 1; i <= width; i++)
            {
                x[IDX(i, j)] = (x0[IDX(i, j)] + a * (x[IDX(i-1, j)] + x[IDX(i+1, j)] + x[IDX(i, j-1)] + x[IDX(i, j+1)]))
                                / (1 + 4 * a);
            }
        }
        setBound(b, x);
    }
}

// Does an advection step on the given components
void FluidSim2D::advect(Component b, float *d, float *d0, float *u, float *v)
{
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dtx, dty;

    dtx = dt * width;
    dty = dt * height;

    for (j = 1; j <= height; j++)
    {
        for (i = 1; i <= width; i++)
        {
            // Apparently this is a "simple linear backtrace"
            x = i - dtx * u[IDX(i, j)];
            y = j - dty * v[IDX(i, j)];
            
            // Keep cells in bounds it looks like
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

            s1 = x - i0;
            s0 = 1 - s1;
            t1 = y - j0;
            t0 = 1 - t1;

            d[IDX(i, j)] = s0 * (t0 * d0[IDX(i0, j0)] + t1 * d0[IDX(i0, j1)]) + 
                           s1 * (t0 * d0[IDX(i1, j0)] + t1 * d0[IDX(i1, j1)]);
        }
    }
    setBound(b, d);
}

// Does some sort of magic projection to get a mass-conserving flow
void FluidSim2D::project(float *u, float *v, float *p, float *div)
{
   int i, j, k;

   float hx, hy;

   hx = hSpacing();
   hy = vSpacing();

   for (j = 1; j <= height; j++)
   {
       for (i = 1; i <= width; i++)
       {
            div[IDX(i, j)] = -0.5 * (hx * u[IDX(i + 1, j)] - hx * u[IDX(i - 1, j)] +
                                     hy * v[IDX(i, j + 1)] - hy * v[IDX(i, j - 1)]);
            p[IDX(i, j)] = 0;
       }
   }
   setBound(NONE, div);
   setBound(NONE, p);

   for (k = 0; k < 20; k++)
   {
       for (j = 1; j <= height; j++)
       {
           for (i = 1; i <= width; i++)
           {
               p[IDX(i, j)] = (div[IDX(i, j)] + p[IDX(i - 1, j)] + p[IDX(i + 1, j)] + p[IDX(i, j - 1)] + p[IDX(i, j + 1)])
                              / 4;
           }
       }
       setBound(NONE, p);
   }

   for (j = 1; j <= height; j++)
   {
       for (i = 1; i <= width; i++)
       {
            u[IDX(i, j)] -= 0.5 * (p[IDX(i + 1, j)] - p[IDX(i - 1, j)]) / hx;
            v[IDX(i, j)] -= 0.5 * (p[IDX(i, j + 1)] - p[IDX(i, j - 1)]) / hy;
       }
   }
   setBound(VERT, u);
   setBound(HORIZ, v);
}

void FluidSim2D::setBound(Component b, float *x)
{
    return;

    int i;

    // Fix edges
    for (i = 1; i <= height; i++)
    {
        x[IDX(0, i)] = (b == VERT) ? -x[IDX(1, i)] : x[IDX(1, i)];
        x[IDX(width + 1, i)] = (b == VERT) ? -x[IDX(width, i)] : x[IDX(width, 1)];
    }

    for (i = 1; i <= width; i++)
    {
        x[IDX(i, 0)] = (b == HORIZ) ? -x[IDX(i, 1)] : x[IDX(i, 1)];
        x[IDX(i, height + 1)] = (b == HORIZ) ? -x[IDX(i, height)] : x[IDX(i, height)];
    }

    // Corner cases
    x[IDX(0, 0)] = 0.5 * (x[IDX(1, 0)] + x[IDX(0, 1)]);
    x[IDX(0, height + 1)] = 0.5 * (x[IDX(1, height + 1)] + x[IDX(0, height)]);
    x[IDX(width + 1, 0)] = 0.5 * (x[IDX(width, 0)] + x[IDX(width + 1, 1)]);
    x[IDX(width + 1, height + 1)] = 0.5 * (x[IDX(width, height + 1)] + x[IDX(width + 1, height)]);
}

// Calculate the size of the grid. Remember that the outside cells are boundary cells.
inline
int FluidSim2D::size()
{
    return (width + 2) * (height + 2);
}

inline
float FluidSim2D::hSpacing()
{
    return 1.0 / (float)width;
}

inline
float FluidSim2D::vSpacing()
{
    return 1.0 / (float)height;
}

void FluidSim2D::debugPrintGrid(float *g)
{
    for (int i = 0; i < width + 2; i++) {
        for (int j = 0; j < height + 2; j++) {
            std::cout << g[IDX(i, j)] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}