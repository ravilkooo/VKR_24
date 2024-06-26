#include "SoundMaterial.h"

///////////////////////////////////////////////
///////////////////////////////////////////////
#define MATERIAL( _material_name, _reflectivity, _scattering ) \
static SoundMaterial new##_material_name##Material()\
{\
	FrequencyBandResponse reflectivity(0), scattering(0);\
	_reflectivity;\
	_scattering;\
	return SoundMaterial( reflectivity, scattering, #_material_name );\
}\

#define NEW_MATERIAL( material_name ) new##material_name##Material()

#define R( f, a ) reflectivity.setFrequency( f, a );
#define S( f, a ) scattering.setFrequency( f, a );

#define CONCRETE_R R( 125.0, 0.99 ) R( 250.0, 0.99 ) R( 500.0, 0.99 ) R( 1000.0, 0.99 ) R( 2000.0, 0.99 ) R( 4000.0, 0.99 )
#define CONCRETE_S S( 125.0, 0.10 ) S( 250.0, 0.11 ) S( 500.0, 0.12 ) S( 1000.0, 0.13 ) S( 2000.0, 0.14 ) S( 4000.0, 0.15 )

MATERIAL(Concrete, CONCRETE_R, CONCRETE_S);
const SoundMaterial SoundMaterial::CONCRETE = NEW_MATERIAL(Concrete);


#define SNOW_R R( 125.0, 0.74 ) R( 250.0, 0.50 ) R( 500.0, 0.32 ) R( 1000.0, 0.22 ) R( 2000.0, 0.22 ) R( 4000.0, 0.22 )
#define SNOW_S S( 125.0, 0.20 ) S( 250.0, 0.30 ) S( 500.0, 0.40 ) S( 1000.0, 0.50 ) S( 2000.0, 0.60 ) S( 4000.0, 0.75 )

MATERIAL(Snow, SNOW_R, SNOW_S);
const SoundMaterial SoundMaterial::SNOW = NEW_MATERIAL(Snow);


#define SILENCE_R R( 125.0, 0.0011 ) R( 250.0, 0.0011 ) R( 500.0, 0.0011 ) R( 1000.0, 0.0011 ) R( 2000.0, 0.0011 ) R( 4000.0, 0.0011 )
#define SILENCE_S S( 125.0, 0.9090909 ) S( 250.0, 0.9090909 ) S( 500.0, 0.9090909 ) S( 1000.0, 0.9090909 ) S( 2000.0, 0.9090909 ) S( 4000.0, 0.9090909 )

MATERIAL(Silence, SILENCE_R, SILENCE_S);
const SoundMaterial SoundMaterial::SILENCE = NEW_MATERIAL(Silence);

///////////////////////////////////////////////
///////////////////////////////////////////////


SoundMaterial::SoundMaterial()
{
	reflectivity = FrequencyBandResponse(0);
	scattering = FrequencyBandResponse(0);
	material_name = "Silence";
	SILENCE_R
	SILENCE_S
}
SoundMaterial::SoundMaterial(FrequencyBandResponse reflectivity, FrequencyBandResponse scattering, std::string material_name):
	reflectivity(reflectivity), scattering(scattering), material_name(material_name)
{
}

std::pair<double, double> SoundMaterial::getClosestParameters(double frequency)
{
	return std::pair<double, double>(reflectivity.getClosestValue(frequency), scattering.getClosestValue(frequency));
}
