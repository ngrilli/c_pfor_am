// Nicolò Grilli
// Università di Bristol
// 10 Febbraio 2024

#include "ACInterfaceAniso.h"

registerMooseObject("c_pfor_amApp", ACInterfaceAniso);

InputParameters
ACInterfaceAniso::validParams()
{
  InputParameters params = ACInterface::validParams();
  params.addClassDescription("Anisotropic gradient energy Allen-Cahn Kernel");
  params.addParam<Real>("e_anisotropy",0.0,"Grain boundary energy anisotropy coefficient. ");
  params.addParam<int>("op",0,"Phase field index used to select the Euler angles. ");
  params.addParam<int>("op_num",0,"Total number of phase fields. ");
  params.addParam<FileName>(
      "Euler_angles_file_name","",
      "Name of the file containing the Euler angles, each row must contain three Euler angles "
      "which correspond to each grain orientation. ");
  return params;
}

ACInterfaceAniso::ACInterfaceAniso(const InputParameters & parameters)
  : ACInterface(parameters),
  _e_anisotropy(getParam<Real>("e_anisotropy")),
  _op(getParam<int>("op")),
  _op_num(getParam<int>("op_num")),
  _Euler_angles_file_name(getParam<FileName>("Euler_angles_file_name")),
  _Euler_angles(0,0,0),
  _R(_Euler_angles)
{
  // Get mobility and kappa derivatives and coupled variable gradients
  for (unsigned int i = 0; i < _n_args; ++i)
  {
    MooseVariable * ivar = _coupled_standard_moose_vars[i];
    const VariableName iname = ivar->name();
    if (iname == _var.name())
    {
      if (isCoupled("args"))
        paramError("args",
                   "The kernel variable should not be specified in the coupled `args` parameter.");
      else
        paramError("coupled_variables",
                   "The kernel variable should not be specified in the coupled `coupled_variables` "
                   "parameter.");
    }

    _dLdarg[i] = &getMaterialPropertyDerivative<Real>("mob_name", i);
    _dkappadarg[i] = &getMaterialPropertyDerivative<Real>("kappa_name", i);
    _d2Ldargdop[i] = &getMaterialPropertyDerivative<Real>("mob_name", iname, _var.name());

    _gradarg[i] = &(ivar->gradSln());

    _d2Ldarg2[i].resize(_n_args);
    for (unsigned int j = 0; j < _n_args; ++j)
      _d2Ldarg2[i][j] = &getMaterialPropertyDerivative<Real>("mob_name", i, j);
  }
  
  if (this->isParamValid("Euler_angles_file_name")) {

    assignEulerAngles();
	  
  } // otherwise Euler angles are zero
  
  computeRotationMatrix();
}

void
ACInterfaceAniso::assignEulerAngles()
{
  // read in the slip system data from auxiliary text file
  MooseUtils::DelimitedFileReader _reader(_Euler_angles_file_name);
  _reader.setFormatFlag(MooseUtils::DelimitedFileReader::FormatFlag::ROWS);
  _reader.read();	
  
  // add check for the size of the input
  if (_reader.getData().size() != _op_num)
    paramError(
        "_op_num",
        "The number of rows in the Euler angles file should match the number of phase fields.");

  for (const auto j : index_range(_reader.getData(_op))) 
  {
    _Euler_angles(j) = _reader.getData(_op)[j];
  }
}

void
ACInterfaceAniso::computeRotationMatrix()
{
  // update function takes angles in degrees
  _R.update(_Euler_angles); // this is passive rotation, see RotationTensor.C
  _crysrot = _R.transpose(); // therefore needs to be transposed
}

Real
ACInterfaceAniso::computeAnisotropy()
{
  std::vector<RealVectorValue> cubic_directions; // in the crystal reference frame
  std::vector<RealVectorValue> rotated_cubic_directions; // in the model reference frame
  Real cos_phi_inc = 0.0; // cosine of inclination angle
  Real temp_abs_cos_phi_inc = 0.0; // and temporary variable to store its absolute value
  Real cos_min_phi_inc = 0.0; // cosine of the minimum angle between grain boundary normal and <001>
  Real sin_min_phi_inc = 1.0; // sine of the minimum angle between grain boundary normal and <001>
  cubic_directions.resize(2*LIBMESH_DIM);
  rotated_cubic_directions.resize(2*LIBMESH_DIM);
  
  // Assign cubic directions and calculate rotated directions
  // checking for only 3 directions may be enough
  for (const auto i : make_range(2*LIBMESH_DIM))
  {
    cubic_directions[i].zero();
    rotated_cubic_directions[i].zero();
  }
  
  for (const auto i : make_range(LIBMESH_DIM))
  {
    cubic_directions[2*i](i) = 1.0;
    cubic_directions[2*i+1](i) = -1.0; 
  }
  
  for (const auto i : make_range(2*LIBMESH_DIM)) {
    for (const auto j : make_range(LIBMESH_DIM)) {
      for (const auto k : make_range(LIBMESH_DIM)) {
        rotated_cubic_directions[i](j) += _crysrot(j,k) * cubic_directions[i](k);
	  }
	}
  }

  // avoid division by zero and check if current grain exists at this _qp
  if (_grad_u[_qp].norm() > 1.0e-12 && _u[_qp] > 1.0e-12) {
	  
    for (const auto i : make_range(2*LIBMESH_DIM)) { // identify the maximum among the cosines
		
      cos_phi_inc = _grad_u[_qp] * rotated_cubic_directions[i] / _grad_u[_qp].norm();
      
      if (std::abs(cos_phi_inc) > temp_abs_cos_phi_inc) { // this <001> direction is closer to GB normal
		  
        temp_abs_cos_phi_inc = std::abs(cos_phi_inc);
        cos_min_phi_inc = cos_phi_inc;
	  }
	}
	
	sin_min_phi_inc = std::sqrt(1.0 - std::pow(cos_min_phi_inc,2));
	  
    return 1.0 + _e_anisotropy * (std::pow(cos_min_phi_inc,4) + std::pow(sin_min_phi_inc,4));
    	  
  } else { // not near a grain boundary, therefore no need for anisotropy 
	  
    return 1.0;
  }
}

Real
ACInterfaceAniso::computeQpResidual()
{
  return computeAnisotropy() * _grad_u[_qp] * kappaNablaLPsi();
}
