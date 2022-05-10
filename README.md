# Binary_Stars_Simulations

I created a python model that simulates orbital motion of hypothetical binary star systems. I created six hypothetical binary star systems by using Newtonian mechanics, and assumed that orbital path is circular for simplicity. I was curious if it was possible to detect a planet orbiting around two stars using radial velocity method which is a method of finding exoplanets by looking at the wobbling motion of the star caused by the gravitational pull of planets orbiting around the star.

First, I had to create a program that calculates the velocity and acceleration of each object and then integrate it to get position and velocity. The calculated values were stored and plotted. I used R to animate some hypothetical systems.

Here are some animations:

![Sim1_noplanet_slow](https://user-images.githubusercontent.com/94130159/167428346-defe9b47-895b-4816-9169-7231ce84774e.gif)

![Sim2_with_planets_slow](https://user-images.githubusercontent.com/94130159/167428386-62119721-f44c-4351-adac-3e331aafd5eb.gif)

![Sim5_noplanet](https://user-images.githubusercontent.com/94130159/167428422-d0608bd7-b1b1-4e7d-8271-d847d2b55d64.gif)


The x-axis velocity of one of the stars in each system was used for spectral analysis. In theory, the gravitatinal pull of the planet would have an effect on star's velocity, allowing us to detect via spectral analysis. However, I was not able to find singificant signals. The mass of the second star is too big and overpowers any other gravitational signals. In order to detect an exoplanet, the exoplanet's mass has to be near the mass of the stars. To detect an exoplanet in binary star systems one might want to use Transit Method.
