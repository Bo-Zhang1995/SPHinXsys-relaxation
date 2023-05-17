/**
 * @file 	relaxation_evolution.cpp
 * @brief 	This is the first case by testing the relaxation with evolution method.
 * @author 	Bo Zhang
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real LL = 1.0;					
Real LH = 1.0;
Real resolution_ref = LH / 50.0;
Real BW = resolution_ref * 2.0;
BoundingBox system_domain_bounds(Vec2d::Zero(), Vec2d(LL, LH));
//----------------------------------------------------------------------
//	Define geometries
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH);
Vec2d water_block_translation = water_block_halfsize;
class Insert : public ComplexShape
{
public:
	explicit Insert(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TransformShape<GeometricShapeBox>>(Transform2d(water_block_halfsize), water_block_translation);
	}
};

int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	sph_system.setRunParticleRelaxation(true);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody body(sph_system, makeShared<Insert>("InsertedBody"));
	/* Change the kernel function to the guass kernel function. */
    //body.sph_adaptation_->resetKernel<KernelTabulated<KernelLaguerreGauss>>(20);
	body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	body.defineParticlesAndMaterial();
	body.addBodyStateForRecording<Vecd>("Position");
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? body.generateParticles<ParticleGeneratorReload>(io_environment, body.getName())
		: body.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation insert_body_inner(body);
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_insert_body_to_vtp(io_environment, { &body });
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(io_environment, { &body });

		/* Relaxation method: including 0th and 1st order consistency. */
		InteractionDynamics<relax_dynamics::CalculateParticleStress> calculate_particle_stress(insert_body_inner, false);
		relax_dynamics::RelaxationStepInner relaxation_0th_inner(insert_body_inner, false);
		relax_dynamics::RelaxationStepImplicitInner relaxation_0th_implicit_inner(insert_body_inner, false);
		relax_dynamics::RelaxationStepByStressInner relaxation_1st_inner(insert_body_inner, false);
		relax_dynamics::RelaxationStepByStressImplicitInner relaxation_1st_implicit_inner(insert_body_inner, false);

		/* Evolution method */
		relax_dynamics::ZeroOrderEvolutionStep zero_order_evolution(insert_body_inner, false);
		relax_dynamics::FirstOrderEvolutionStep first_order_evolution(insert_body_inner, false);
		InteractionSplit<relax_dynamics::CorrectionMatrixRegularization> correction_matrix_regularization(insert_body_inner);

		PeriodicConditionUsingCellLinkedList periodic_condition_x(body, body.getBodyShapeBounds(), xAxis);
		PeriodicConditionUsingCellLinkedList periodic_condition_y(body, body.getBodyShapeBounds(), yAxis);

		/* Update relaxation residue. */
		InteractionDynamics<relax_dynamics::CheckCorrectedZeroOrderConsistency> check_corrected_zero_order_consistency(insert_body_inner);
		ReduceAverage<QuantitySummation<Real>> calculate_particle_average_zero_error(body, "corrected_zero_order_error");
		ReduceDynamics<QuantityMaximum<Real>> calculate_particle_maximum_zero_error(body, "corrected_zero_order_error");
		InteractionDynamics<relax_dynamics::CheckCorrectedFirstOrderConsistency> check_corrected_first_order_consistency(insert_body_inner);
		ReduceAverage<QuantitySummation<Real>> calculate_particle_average_first_error(body, "corrected_first_order_error");
		ReduceDynamics<QuantityMaximum<Real>> calculate_particle_maximum_first_error(body, "corrected_first_order_error");
		//----------------------------------------------------------------------  
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_insert_body_particles.exec(0.25);
		sph_system.initializeSystemCellLinkedLists();
		periodic_condition_x.update_cell_linked_list_.exec();
		periodic_condition_y.update_cell_linked_list_.exec();
		sph_system.initializeSystemConfigurations();
		write_insert_body_to_vtp.writeToFile(0);

		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		TickCount t1 = TickCount::now();

		int ite = 0; //iteration step for the total relaxation step.
		int ite_p = 0; //iteration step for the 0th order relaxation step.
		int ite_s = 0; //iteration step for the 1st order relaxation step.
		int ite_loop = 0; //iteration loop for the whole relaxation process.

		Real last_zero_maximum_residual = 1;
		Real last_zero_average_residual = 1;
		Real last_first_maximum_residual = 1;
		Real last_first_average_residual = 1;

		Real current_zero_maximum_residual = 1; //maximum zero order consistency residual.
		Real current_zero_average_residual = 1; //average zero order consistency residual.
		Real current_first_maximum_residual = 1; //maximum first order consistency residual.
		Real current_first_average_residual = 1; //average first order consistency residual.

		GlobalStaticVariables::physical_time_ = ite;

		std::string average_zero_error = io_environment.output_folder_ + "/" + "average_zero_error.dat";
		std::ofstream out_average_zero_error(average_zero_error.c_str(), std::ios::app);

		std::string maximum_zero_error = io_environment.output_folder_ + "/" + "maximum_zero_error.dat";
		std::ofstream out_maximum_zero_error(maximum_zero_error.c_str(), std::ios::app);

		std::string average_first_error = io_environment.output_folder_ + "/" + "average_first_error.dat";
		std::ofstream out_average_first_error(average_first_error.c_str(), std::ios::app);

		std::string maximum_first_error = io_environment.output_folder_ + "/" + "maximum_first_error.dat";
		std::ofstream out_maximum_first_error(maximum_first_error.c_str(), std::ios::app);

		std::string evolution_zero_step = io_environment.output_folder_ + "/" + "evolution_zero_step.dat";
		std::ofstream out_evolution_zero_step(evolution_zero_step.c_str(), std::ios::app);

		std::string evolution_first_step = io_environment.output_folder_ + "/" + "evolution_first_step.dat";
		std::ofstream out_evolution_first_step(evolution_first_step.c_str(), std::ios::app);

		/* The procedure to obtain uniform particle distribution that satisfies the 0th order consistency. */
		while (current_zero_maximum_residual > 0.0001)
		//while (current_zero_average_residual > 0.001 || current_zero_maximum_residual > 0.005)
		{
			periodic_condition_x.bounding_.exec();
			periodic_condition_y.bounding_.exec();
			body.updateCellLinkedList();
			periodic_condition_x.update_cell_linked_list_.exec();
			periodic_condition_y.update_cell_linked_list_.exec();
			insert_body_inner.updateConfiguration();

			// Zero order evolution procedure.
			//relaxation_0th_implicit_inner.exec(0.1);
			// First + zero order evolution procedure.
			calculate_particle_stress.exec();
			relaxation_1st_implicit_inner.exec(0.1);

			ite++;

			if (ite % 100 == 0)
			{
				check_corrected_zero_order_consistency.exec();
				current_zero_average_residual = calculate_particle_average_zero_error.exec();
				current_zero_maximum_residual = calculate_particle_maximum_zero_error.exec();

				std::cout << std::fixed << std::setprecision(9) << "0th relaxation steps for the body N = " << ite << "\n";
				std::cout << "$$ The 0th consistency error: maximum = " << current_zero_maximum_residual << "; average = " << current_zero_average_residual << std::endl;
				write_insert_body_to_vtp.writeToFile(ite);
			}
		}
		/*calculate_particle_stress.exec();

		check_corrected_first_order_consistency.exec();
		current_first_average_residual = calculate_particle_average_first_error.exec();
		current_first_maximum_residual = calculate_particle_maximum_first_error.exec();
		std::cout << "¡ïThe 1st consistency error after 0th relaxation: maximum = " << current_first_maximum_residual << "; average = " << current_first_average_residual << std::endl;
		*/

		check_corrected_zero_order_consistency.exec();
		current_zero_average_residual = calculate_particle_average_zero_error.exec();
		current_zero_maximum_residual = calculate_particle_maximum_zero_error.exec();
		std::cout << "@The 0th consistency error: maximum = " << current_zero_maximum_residual << "; average = " << current_zero_average_residual << std::endl;

		ite++;
		write_insert_body_to_vtp.writeToFile(ite);

		///* The procedure to test correction matrix evolution that satisfies the 0th order consistency. */
		//while (current_first_average_residual > 0.005)
		////while (current_zero_average_residual > 0.001 || current_zero_maximum_residual > 0.005)
		//{
		//	periodic_condition_x.bounding_.exec();
		//	periodic_condition_y.bounding_.exec();
		//	body.updateCellLinkedList();
		//	periodic_condition_x.update_cell_linked_list_.exec();
		//	periodic_condition_y.update_cell_linked_list_.exec();
		//	insert_body_inner.updateConfiguration();

		//	// first order evolution procedure. 
		//	first_order_evolution.exec(0.1);
		//	ite++;

		//	if (ite % 100 == 0)
		//	{
		//		//correction_matrix_regularization.exec(0.1);
		//		check_corrected_first_order_consistency.exec();
		//		current_first_average_residual = calculate_particle_average_first_error.exec();
		//		current_first_maximum_residual = calculate_particle_maximum_first_error.exec();

		//		std::cout << std::fixed << std::setprecision(9) << "1st relaxation steps for the body N = " << ite << "\n";
		//		std::cout << "$$ The 1st consistency error: maximum = " << current_first_maximum_residual << "; average = " << current_first_average_residual << std::endl;
		//		write_insert_body_to_vtp.writeToFile(ite);
		//	}
		//}

		///* The initial error for the 0th and 1st order consistency.*/
		//check_corrected_zero_order_consistency.exec();
		//last_zero_maximum_residual = calculate_particle_maximum_zero_error.exec();
		//last_zero_average_residual = calculate_particle_average_zero_error.exec();
		//current_zero_maximum_residual = last_zero_maximum_residual;
		//current_zero_average_residual = last_zero_average_residual;

		//check_corrected_first_order_consistency.exec();
		//last_first_maximum_residual = calculate_particle_maximum_first_error.exec();
		//last_first_average_residual = calculate_particle_average_first_error.exec();
		//current_first_maximum_residual = last_first_maximum_residual;
		//current_first_average_residual = last_first_average_residual;

		//out_average_zero_error << std::fixed << std::setprecision(12) << ite << "    " << current_zero_average_residual << "\n";
		//out_maximum_zero_error << std::fixed << std::setprecision(12) << ite << "    " << current_zero_maximum_residual << "\n";
		//out_average_first_error << std::fixed << std::setprecision(12) << ite << "    " << current_first_average_residual << "\n";
		//out_maximum_first_error << std::fixed << std::setprecision(12) << ite << "    " << current_first_maximum_residual << "\n";

		///* The termination condition is the maximum and average 0th and 1st consistency errors are all smaller than 0.001. */
		////while ((current_zero_average_residual > 0.001 || current_zero_maximum_residual > 0.005 || current_first_average_residual > 0.001 || current_first_maximum_residual > 0.005) && (ite_loop < 1000))
		//while ((current_zero_average_residual > 0.001 || current_first_average_residual > 0.007) && (ite_loop < 1000))
		//{
		//	std::cout << "******************" << ite_loop << " loop start******************" << std::endl;

		//	/*****************************************************************************/
		//	/* Evolution of correction matrix doesn't change the particle configuration. */
		//	periodic_condition_x.bounding_.exec();
		//	periodic_condition_y.bounding_.exec();
		//	body.updateCellLinkedList();
		//	periodic_condition_x.update_cell_linked_list_.exec();
		//	periodic_condition_y.update_cell_linked_list_.exec();
		//	insert_body_inner.updateConfiguration();

		//	//while ((current_first_maximum_residual > 0.99 * last_first_maximum_residual) || (current_first_average_residual > 0.99 * last_first_average_residual))
		//	while((current_first_average_residual > 0.99 * last_first_average_residual) && (current_first_average_residual > 0.007))
		//	{
		//		first_order_evolution.exec(0.0001);
		//		ite++; ite_s++;

		//		if (ite % 100 == 0)
		//		{
		//			/* regularize the correction matrix */
		//			correction_matrix_regularization.exec(0.0001);

		//			check_corrected_first_order_consistency.exec();
		//			current_first_average_residual = calculate_particle_average_first_error.exec();
		//			current_first_maximum_residual = calculate_particle_maximum_first_error.exec();
		//			out_average_first_error << std::fixed << std::setprecision(12) << ite << "    " << current_first_average_residual << "\n";
		//			out_maximum_first_error << std::fixed << std::setprecision(12) << ite << "    " << current_first_maximum_residual << "\n";

		//			std::cout << std::fixed << std::setprecision(9) << "¡ïThis is 1st relaxation steps for the body N = " << ite << "\n";
		//			std::cout << "¡ïThe 1st consistency error: maximum = " << current_first_maximum_residual << "; average = " << current_first_average_residual << std::endl;
		//			write_insert_body_to_vtp.writeToFile(ite);
		//		}
		//	}
		//	out_evolution_first_step << std::fixed << std::setprecision(12) << ite << "    " << ite_s << "\n";
		//	//std::cout << "¡ïThis loop for 1st relaxation is finished: maximum = " << current_first_maximum_residual << "; average = " << current_first_average_residual << std::endl;

		//	ite_s = 0; //reset the number of iteration of 1st consistency.
		//	last_first_average_residual = current_first_average_residual;
		//	last_first_maximum_residual = current_first_maximum_residual;
		//	write_insert_body_to_vtp.writeToFile(ite);
		//	
		//	/* Check the 0th order consistency error after 1st evolution. */
		//	check_corrected_zero_order_consistency.exec();
		//	current_zero_average_residual = calculate_particle_average_zero_error.exec();
		//	current_zero_maximum_residual = calculate_particle_maximum_zero_error.exec();
		//	out_average_zero_error << std::fixed << std::setprecision(12) << ite << "    " << current_zero_average_residual << "\n";
		//	out_maximum_zero_error << std::fixed << std::setprecision(12) << ite << "    " << current_zero_maximum_residual << "\n";
		//	std::cout << "@The 0th consistency error after 1st relaxation: maximum = " << current_zero_maximum_residual << "; average = " << current_zero_average_residual << std::endl;

		//	/*****************************************************************/
		//	/* Evolution of position does change the particle configuration. */
		//	//while ((current_zero_average_residual > 0.99 * last_zero_average_residual) || (current_zero_maximum_residual > 0.99 * last_zero_maximum_residual))
		//	//while(current_zero_average_residual > 0.001 || current_zero_maximum_residual > 0.005)
		//	//while (current_zero_average_residual > 0.005)
		//	while((current_zero_average_residual > 0.99 * last_zero_average_residual) && (current_zero_average_residual > 0.001))
		//	{
		//		periodic_condition_x.bounding_.exec();
		//		periodic_condition_y.bounding_.exec();
		//		body.updateCellLinkedList();
		//		periodic_condition_x.update_cell_linked_list_.exec();
		//		periodic_condition_y.update_cell_linked_list_.exec();
		//		insert_body_inner.updateConfiguration();

		//		/* Zero order evolution procedure. */
		//		zero_order_evolution.exec(0.5);
		//		ite++; ite_p++;

		//		if (ite % 100 == 0)
		//		{
		//			check_corrected_zero_order_consistency.exec();
		//			current_zero_average_residual = calculate_particle_average_zero_error.exec();
		//			current_zero_maximum_residual = calculate_particle_maximum_zero_error.exec();
		//			out_average_zero_error << std::fixed << std::setprecision(12) << ite << "   " << current_zero_average_residual << "\n";
		//			out_maximum_zero_error << std::fixed << std::setprecision(12) << ite << "   " << current_zero_maximum_residual << "\n";

		//			std::cout << std::fixed << std::setprecision(9) << "@This is the 0th relaxation steps for the body N = " << ite << "\n";
		//			std::cout << "@The 0th consistency error: maximum = " << current_zero_maximum_residual << "; average = " << current_zero_average_residual << std::endl;
		//			write_insert_body_to_vtp.writeToFile(ite);
		//		}
		//	}
		//	out_evolution_zero_step << std::fixed << std::setprecision(12) << ite << "    " << ite_p << "\n";
		//	std::cout << "@This loop for 0th relaxation is finished: maximum = " << current_zero_maximum_residual << "; average = " << current_zero_average_residual << std::endl;

		//	ite_p = 0; //reset the number of iteration of 0th consistency.
		//	last_zero_average_residual = current_zero_average_residual;
		//	last_zero_maximum_residual = current_zero_maximum_residual;
		//	write_insert_body_to_vtp.writeToFile(ite);

		//	///* Final check the 0th consistency error. */
		//	//check_corrected_zero_order_consistency.exec();
		//	//current_zero_average_residual = calculate_particle_average_zero_error.exec();
		//	//current_zero_maximum_residual = calculate_particle_maximum_zero_error.exec();
		//	//out_average_zero_error << std::fixed << std::setprecision(12) << ite << "    " << current_zero_average_residual << "\n";
		//	//out_maximum_zero_error << std::fixed << std::setprecision(12) << ite << "    " << current_zero_maximum_residual << "\n";
		//	//std::cout << "@@ The 0th consistency error after 1st relaxation: maximum = " << current_zero_maximum_residual << "; average = " << current_zero_average_residual << std::endl;

		//	/* Final check the 1st consistency error. */
		//	check_corrected_first_order_consistency.exec();
		//	current_first_average_residual = calculate_particle_average_first_error.exec();
		//	current_first_maximum_residual = calculate_particle_maximum_first_error.exec();
		//	out_average_first_error << std::fixed << std::setprecision(12) << ite << "    " << current_first_average_residual << "\n";
		//	out_maximum_first_error << std::fixed << std::setprecision(12) << ite << "    " << current_first_maximum_residual << "\n";
		//	std::cout << "¡ïThe 1st consistency error after 0th relaxation: maximum = " << current_first_maximum_residual << "; average = " << current_first_average_residual << std::endl;

		//	ite_loop++;
		//}

		std::cout << "The final 0th error: maximum = " << current_zero_maximum_residual << "; average = " << current_zero_average_residual << std::endl;
		std::cout << "The final 1st error: maximum = " << current_first_maximum_residual << "; average = " << current_first_average_residual << std::endl;
		std::cout << "The physical relaxation process of body finish !" << std::endl;
		
		/* Output results. */
		write_particle_reload_files.writeToFile(0);
		TickCount t2 = TickCount::now();
		TickCount::interval_t tt;
		tt = t2 - t1;
		std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
		return 0;
	}
}
