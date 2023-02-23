#include <deltaX.h> 
#include <aspect/boundary_temperature/constant.h>
#include <aspect/lateral_averaging.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      deltaX<dim>::
      deltaX ()
        :
        DataPostprocessorScalar<dim> ("deltaX", // to be changed
                                      update_values | update_quadrature_points )
      {
      }
      
      template <int dim>
      void
      deltaX<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const double max_depth = this->get_geometry_model().maximal_depth();
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());*

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const double temperature = input_data.solution_values[q][this->introspection().component_indices.temperature];
            const double depth = this->get_geometry_model().depth (input_data.evaluation_points[q]);
            const double gravity = 1.4;
            const double density = 916;
            const double pressure  = gravity * density * depth;
            const double T_liquidus = 73.2 * pow(1 - pressure/395,1./9.);

            const double deltaX = 1500 * (temperature - T_liquidus); 
            computed_quantities[q](0) = deltaX;
            }
      }

      template <int dim>
      void
      TemperatureAnomaly<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("deltaX");
            {
              prm.declare_entry ("Number of depth slices","20",
                                 Patterns::Integer (1),
                                 "Number of depth slices used to define "
                                 "average temperature.");
              prm.declare_entry ("Use maximal temperature for bottom","true",
                                 Patterns::Bool(),
                                 "If true, use the specified boundary temperatures as average temperatures at the surface. "
                                 "If false, extrapolate the temperature gradient between the first and second cells to the surface. "
                                 "This option will only work for models with a fixed surface temperature. ");
              prm.declare_entry ("Use minimal temperature for surface","true",
                                 Patterns::Bool(),
                                 "Whether to use the minimal specified boundary temperature as the bottom boundary temperature. "
                                 "This option will only work for models with a fixed bottom boundary temperature. ");

            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      template <int dim>
      void
      TemperatureAnomaly<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("deltaX");
            {
              n_slices = prm.get_integer("Number of depth slices");
              extrapolate_surface = !prm.get_bool("Use minimal temperature for surface");
              extrapolate_bottom = !prm.get_bool("Use maximal temperature for bottom");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(deltaX,
                                                  "deltaX",
                                                  "A visualization output postprocessor that outputs the temperature minus the depth-average of the temperature."
                                                  "The average temperature is calculated using the lateral averaging function from the ``depth average'' "
                                                  "postprocessor and interpolated linearly between the layers specified through ``Number of depth slices''")
    }
  }
}

