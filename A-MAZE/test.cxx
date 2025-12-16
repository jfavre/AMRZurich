#include <map>
#include <string>
int main()
{
  std::map<std::string, std::string> PVlabels;

  PVlabels["Density"] = "Log10(Density) [gr/cm^3]";
  PVlabels["Velocity"] = "Velocity [cm/s]";
  PVlabels["Energy Density"] = "Energy Density [erg/cm^3]";
  PVlabels["Thermal Energy Density"] = "Thermal Energy Density [erg/cm^3]";
  PVlabels["Pressure"] = "Log10(Pressure) [erg/cm^3]";
  PVlabels["Mach Number"] = "Mach Number";
  PVlabels["Temperature"] = "Log10(Temperature) [Kelvin]";

  PVlabels.clear();
}
