import matplotlib.pyplot as plt
from rockets import Rocket

class FlightCalc:
    def plot_tank_height_to_attitude(self, tank_heights=[1, 2, 3],
                                     tank_rads=[i for i in range(1, 15)]):
        mdot = Rocket().mdot
        rocket_heights = []
        last_rads, masses, prop_masses, max_heights = [], [], [], []
        for h in tank_heights:
            rocs = []
            for r in tank_rads:
                roc = Rocket(radius=r / 100, tanks_height=h)
                if roc.check():
                    rocs.append(roc)
                else:
                    break
            flights = [roc.max_height() for roc in rocs]
            rads = [roc.radius for roc in rocs]
            prop_masses.append(round(rocs[-1].tanks.prop_mass, 2))
            masses.append(round(rocs[-1].total_mass(), 2))
            rocket_heights.append(round(rocs[-1].base_h, 2))
            last_rads.append(rocs[-1].radius)
            max_heights.append(round(flights[-1], 2))

            plt.plot(rads, flights, label=f"Tanks height: {h} m.", marker="o")

        info_dict = {
            "attitudes": max_heights,
            "rocket heights": rocket_heights,
            "tank heights": tank_heights,
            "rocket radius": last_rads,
            "masses": masses,
            "propellant masses": prop_masses,
                }
        print(info_dict)

        plt.figtext(0.1, 0.8, f"Mass Flow {mdot}")
        plt.title("Max heights for Flight with Drag of the Rocket.")
        plt.legend()
        plt.xlabel("Rocket radius, cm")
        plt.ylabel("Max height, meters")
        plt.show()
        plt.close()
        return info_dict


calc = FlightCalc().plot_tank_height_to_attitude()