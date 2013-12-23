#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include "FluidSim2D.h"
#include "FluidSim3D.h"
#include <glext.h>

#define FLUID3D

int main()
{
    int winWidth = 500;
    int winHeight = 500;

    // create the window
    sf::Window window(sf::VideoMode(winWidth, winHeight), "Fluid Simulation", sf::Style::Default, sf::ContextSettings(32));
    window.setVerticalSyncEnabled(true);

    // Limit to 60 fps
    window.setFramerateLimit(60);

    // Fluid sim init
    
#ifndef FLUID3D
    FluidSim2D simulator(192 / 2, 108 / 2, 0.1, 0.0, 0.0);
#else
    FluidSim3D simulator(24, 24, 24, 0.1, 0.0, 0.0, false);
#endif

    // OpenGL initialization
    glClearColor(0, 0, 0, 0);
    glEnable(GL_TEXTURE_3D);

    // run the main loop
    bool running = true;
    while (running)
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // handle events
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
            {
                // end the program
                running = false;
            }
            else if (event.type == sf::Event::Resized)
            {
                // adjust the viewport when the window is resized
            }
            else if (event.type == sf::Event::MouseWheelMoved)
            {
#ifdef FLUID3D
                simulator.dist += event.mouseWheel.delta * -0.1;
#endif
            }
        }

        simulator.reinit();

        // Collect UI events
        if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
            simulator.setSource(sf::Mouse::getPosition(window), winWidth, winHeight);
        if (sf::Mouse::isButtonPressed(sf::Mouse::Right))
            simulator.setVelocity(sf::Mouse::getPosition(window), winWidth, winHeight, 5.0f);

#ifdef FLUID3D
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::A))
            simulator.yrot += 1;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::D))
            simulator.yrot -= 1;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::W))
            simulator.xrot += 1;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::S))
            simulator.xrot -= 1;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::B))
            simulator.bounds = !simulator.bounds;
#endif

        simulator.simulate();

        simulator.draw(winWidth, winHeight);

        // end the current frame (internally swaps the front and back buffers)
        window.display();
    }

    // release resources...

    return 0;
}