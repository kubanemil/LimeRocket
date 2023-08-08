# LimeRocket. Making Optimal Liquid-Propellant Rockets

## How Application Works

[LimeRocket](https://github.com/kubanemil/rocket) is a Python **framework** that uses modern **rocket science** theory to model the liquid propellant rocket. This includes:

-   Finding the optimal configuration and sizing of the rocket
-   Modeling the rocket’s engine.
-   Calculation of variables that used to calculate the **efficiency** of the flight.

LimeRocket uses [NASA’s CEA](https://cearun.grc.nasa.gov/) (Chemical Equilibrium Application) to calculate and model the rocket’s engine.

Flight simulations based on proven dynamic equations.

![](media/602bf7cc4894d10f76fd83b8cff4e0e0.gif)

![](media/dfb486d2d7f4fdb1bfacb63f3196d34c.gif)

Also, LimeRocket finds optimal rocket version for maximum height fly:

![](media/2a081690ce6e8268981b8c6ea93ed4f5.png)

## Why?

For enthusiastic people, who want to build their own advanced and modern rocket, there is yet no tool that can help them with complex calculations. Although there are some tools that provide assistance, like [OpenRocket](https://openrocket.info/) (for solid-propellant, model rockets) and [NASA’s CEARUN](https://cearun.grc.nasa.gov/) (for modeling the engine), their function are narrow and can’t serve as a framework that provides broad opportunities for building LP (liquid-propellant) rockets.

LimeRocket serves exactly that purpose: it can simultaneously model the flight and configuration of the LP rocket. That is why it is perfect for those people, who begin to build their own rockets!

**How to Use it**
1. Download the code: 
``` git clone https://github.com/kubanemil/LimeRocket.git ```
2. Install all of the requirements by: 
``` pip install -r requirements.txt ```
3. Get to rockets.py file and create the Rocket() class instance:
``` roc = Rocket() ```
4. Then use .plot_rocket() method:
``` roc.plot_rocket() ```
5. To see the engine type:
```roc.plot_engine()```

## Variable Explanation

### Rocket() class
 - <b>tanks_height</b> - combined height of oxidizer tank and fuel tank (in meters)
 - <b>o_f</b> - mass ratio of oxidizer to fuel (ox_mass/fuel_mass)
 - <b>radius</b> = outer radios of rocket's base (in meters)
 - <b>ho2_per</b> = concentration of HO2 in percents (0% - 100%)
 - <b>others_mass</b> - mass of the scheme and other minor details.
 - <b>KPD</b> - relative efficiency of engine in ratio (0-1.0)
 - <b>mdot</b> - the combined mass flow of fuel+oxidizer (kg/sec)
 - <b>ch_pressure</b> - pressure inside the engine's chamber.
 - <b>Po_Pf</b> - ratio of initial pressure inside the tanks to the final.
 - <b>Cd</b> - drag coefficient of the Rocket.

### Tank() class
 - <b>radius</b> - 
 - <b></b> - 
 - <b></b> - 
 - <b></b> - 
 - <b></b> - 
 - <b></b> - 
 - <b>OOO</b> - 
 - 