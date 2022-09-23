#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 2.0;						   /**< Tank length. */
Real DH = 10.0;						   /**< Tank height. */
Real particle_spacing_ref = DL / 40.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;	   /**< Extending width for BCs. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f1 = 1.0;								  /**< Reference density of the bottom drop of water. */
Real rho0_f2 = 1.0;							  /**< Reference density of the upper drop of water. */
Real gravity_1 = 1.0;							  /**< Gravity force of the bottom drop of water. */
Real gravity_2 = -1.0;							  /**< Gravity force of the upper drop of water. */
Real U_max = 1.0;								  /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;						  /**< Reference sound speed. */
Real mu_f1 = 5.0e-2;								  /**< Lower drop viscosity. */
Real mu_f2 = 5.0e-2;								  /**< Upper drop viscosity. */
Real contact_angle = (150.0 / 180.0) * 3.1415926; /**< Contact angle with Wall. */
Real tension_force = 0.008;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
/** create two water block shapes */

std::vector<Vecd> createLowerWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.375 * DL, 0.0));
	water_block_shape.push_back(Vecd(0.375 * DL, 0.35));
	water_block_shape.push_back(Vecd(0.625 * DL, 0.35));
	water_block_shape.push_back(Vecd(0.625 * DL, 0.0));
	water_block_shape.push_back(Vecd(0.375 * DL, 0.0));
	return water_block_shape;
}
std::vector<Vecd> createUpperWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.375 * DL, DH-0.35));
	water_block_shape.push_back(Vecd(0.375 * DL, DH));
	water_block_shape.push_back(Vecd(0.625 * DL, DH));
	water_block_shape.push_back(Vecd(0.625 * DL, DH-0.35));
	water_block_shape.push_back(Vecd(0.375 * DL, DH-0.35));
	return water_block_shape;
}
/** create outer wall shape */
std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, -BW));

	return outer_wall_shape;
}
/** create inner wall shape */
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(0.0, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, DH));
	inner_wall_shape.push_back(Vecd(DL, DH));
	inner_wall_shape.push_back(Vecd(DL, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, 0.0));

	return inner_wall_shape;
}
//----------------------------------------------------------------------
//	Water block body with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
class LowerWaterBlock : public FluidBody
{
public:
	LowerWaterBlock(SPHSystem &sph_system, const string &body_name)
		: FluidBody(sph_system, body_name, makeShared<SPHAdaptation>(1.3, 1))
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createLowerWaterBlockShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};

class UpperWaterBlock : public FluidBody
{
public:
	UpperWaterBlock(SPHSystem &sph_system, const string &body_name)
		: FluidBody(sph_system, body_name, makeShared<SPHAdaptation>(1.3, 1))
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createUpperWaterBlockShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};

//----------------------------------------------------------------------
//	Wall boundary body definition.
//----------------------------------------------------------------------
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, const std::string &body_name)
		: SolidBody(sph_system, body_name, makeShared<SPHAdaptation>(1.3, 1))
	{
		/** Geomtry definition. */
		std::vector<Vecd> outer_shape = createOuterWallShape();
		std::vector<Vecd> inner_shape = createInnerWallShape();
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(outer_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_shape, ShapeBooleanOps::sub);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
