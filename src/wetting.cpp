#include "wetting.h"
#include "sphinxsys.h"
using namespace SPH;

int main()
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	sph_system.restart_step_ = 0;
	/** I/O environment. */
	In_Output in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	LowerWaterBlock water_block1(sph_system, "WaterBody1");
	FluidParticles water_particles1(water_block1, makeShared<WeaklyCompressibleFluid>(rho0_f1, c_f, mu_f1));

	UpperWaterBlock water_block2(sph_system, "WaterBody2");
	FluidParticles water_particles2(water_block2, makeShared<WeaklyCompressibleFluid>(rho0_f2, c_f, mu_f2));

	WallBoundary wall_boundary(sph_system, "Wall");
	SolidParticles wall_particles(wall_boundary);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexBodyRelation water_water_complex1(water_block1, {&water_block2});
	BodyRelationContact water_wall_contact1(water_block1, {&wall_boundary});
	ComplexBodyRelation water_water_complex2(water_block2, {&water_block1});
	BodyRelationContact water_wall_contact2(water_block2, {&wall_boundary});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Define external force. */
	Gravity gravity1(Vecd(0.0, gravity_1));
	Gravity gravity2(Vecd(0.0, gravity_2));
	/** Initialize particle acceleration. */
	TimeStepInitialization initialize_a_water_step1(water_block1, gravity1);
	TimeStepInitialization initialize_a_water_step2(water_block2, gravity2);
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationFreeSurfaceComplex
		update_water_density_by_summation1(water_water_complex1.inner_relation_, water_wall_contact1);
	fluid_dynamics::DensitySummationFreeSurfaceComplex
		update_water_density_by_summation2(water_water_complex2.inner_relation_, water_wall_contact2);
	/** Time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_water_advection_time_step_size1(water_block1, U_max);
	fluid_dynamics::AdvectionTimeStepSize get_water_advection_time_step_size2(water_block2, U_max);
	/** Time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_water_time_step_size1(water_block1);
	fluid_dynamics::AcousticTimeStepSize get_water_time_step_size2(water_block2);
	/** Pressure relaxation for water by using position verlet time stepping. */
	fluid_dynamics::PressureRelaxationRiemannWithWall
		water_pressure_relaxation1(water_water_complex1.inner_relation_, water_wall_contact1);
	fluid_dynamics::PressureRelaxationRiemannWithWall
		water_pressure_relaxation2(water_water_complex2.inner_relation_, water_wall_contact2);
	fluid_dynamics::DensityRelaxationRiemannWithWall
		water_density_relaxation1(water_water_complex1.inner_relation_, water_wall_contact1);
	fluid_dynamics::DensityRelaxationRiemannWithWall
		water_density_relaxation2(water_water_complex2.inner_relation_, water_wall_contact2);
	/** Viscous acceleration. */
	fluid_dynamics::ViscousAccelerationMultiPhase
		water_viscous_acceleration1(water_water_complex1);
	fluid_dynamics::ViscousAccelerationMultiPhase
		water_viscous_acceleration2(water_water_complex2);
	/** Suface tension and wetting effects. */
	fluid_dynamics::FreeSurfaceIndicationComplex
		surface_detection1(water_water_complex1.inner_relation_, water_wall_contact1);
	fluid_dynamics::FreeSurfaceIndicationComplex
		surface_detection2(water_water_complex2.inner_relation_, water_wall_contact2);
	fluid_dynamics::ColorFunctionGradientComplex
		color_gradient1(water_water_complex1.inner_relation_, water_wall_contact1);
	fluid_dynamics::ColorFunctionGradientComplex
		color_gradient2(water_water_complex2.inner_relation_, water_wall_contact2);
	fluid_dynamics::ColorFunctionGradientInterplationInner
		color_gradient_interpolation1(water_water_complex1.inner_relation_);
	fluid_dynamics::ColorFunctionGradientInterplationInner
		color_gradient_interpolation2(water_water_complex2.inner_relation_);
	fluid_dynamics::SurfaceTensionAccelerationInner
		surface_tension_acceleration1(water_water_complex1.inner_relation_, tension_force);
	fluid_dynamics::SurfaceTensionAccelerationInner
		surface_tension_acceleration2(water_water_complex2.inner_relation_, tension_force);
	/** Wetting effects. */
	fluid_dynamics::SurfaceNormWithWall
		wetting_norm1(water_wall_contact1, contact_angle);
	fluid_dynamics::SurfaceNormWithWall
		wetting_norm2(water_wall_contact2, contact_angle);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	/** Output the body states. */
	BodyStatesRecordingToVtp body_states_recording(in_output, sph_system.real_bodies_);
	/** Output the body states for restart simulation. */
	RestartIO restart_io(in_output, sph_system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromBodyShape();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	/** If the starting time is not zero, please setup the restart time step ro read in restart states. */
	if (sph_system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.restart_step_);
		water_block1.updateCellLinkedList();
		water_block2.updateCellLinkedList();
		water_water_complex1.updateConfiguration();
		water_wall_contact1.updateConfiguration();
		water_water_complex2.updateConfiguration();
		water_wall_contact2.updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = sph_system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 30.0;		 /**< End time. */
	Real D_Time = End_Time / 500; /**< Time stamps for output of body states. */
	Real dt = 0.0;				 /**< Default acoustic time step sizes. */
	/** statistics for computing CPU time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	tick_count::interval_t interval_computing_time_step;
	tick_count::interval_t interval_computing_pressure_relaxation;
	tick_count::interval_t interval_updating_configuration;
	tick_count time_instance;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	body_states_recording.writeToFile();
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** Acceleration due to viscous force and gravity. */
			time_instance = tick_count::now();
			initialize_a_water_step1.parallel_exec();
			initialize_a_water_step2.parallel_exec();

			Real Dt_f1 = get_water_advection_time_step_size1.parallel_exec();
			Real Dt_f2 = get_water_advection_time_step_size2.parallel_exec();
			Real Dt = SMIN(Dt_f1, Dt_f2);

			update_water_density_by_summation1.parallel_exec();
			update_water_density_by_summation2.parallel_exec();

			water_viscous_acceleration1.parallel_exec();
			water_viscous_acceleration2.parallel_exec();

			surface_detection1.parallel_exec();
			surface_detection2.parallel_exec();
			color_gradient1.parallel_exec();
			color_gradient2.parallel_exec();
			color_gradient_interpolation1.parallel_exec();
			color_gradient_interpolation2.parallel_exec();
			wetting_norm1.parallel_exec();
			wetting_norm2.parallel_exec();
			surface_tension_acceleration1.parallel_exec();
			surface_tension_acceleration2.parallel_exec();

			interval_computing_time_step += tick_count::now() - time_instance;

			/** Dynamics including pressure relaxation. */
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				Real dt_f1 = get_water_time_step_size1.parallel_exec();
				Real dt_f2 = get_water_time_step_size2.parallel_exec();
				dt = SMIN(SMIN(dt_f2, dt_f2), Dt);

				water_pressure_relaxation1.parallel_exec(dt);
				water_pressure_relaxation2.parallel_exec(dt);

				water_density_relaxation1.parallel_exec(dt);
				water_density_relaxation2.parallel_exec(dt);

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			interval_computing_pressure_relaxation += tick_count::now() - time_instance;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			/** Update cell linked list and configuration. */
			time_instance = tick_count::now();

			water_block1.updateCellLinkedList();
			water_water_complex1.updateConfiguration();
			water_wall_contact1.updateConfiguration();

			water_block2.updateCellLinkedList();
			water_water_complex2.updateConfiguration();
			water_wall_contact2.updateConfiguration();

			interval_updating_configuration += tick_count::now() - time_instance;
		}

		tick_count t2 = tick_count::now();
		body_states_recording.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}

	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
			  << " seconds." << std::endl;
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
			  << interval_computing_time_step.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_computing_pressure_relaxation = "
			  << interval_computing_pressure_relaxation.seconds() << "\n";
	std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
			  << interval_updating_configuration.seconds() << "\n";

	return 0;
}
