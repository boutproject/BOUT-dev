// Class to deal with sources and source parameters
class Sources {
public:
  Sources();
  ~Sources();
  void initialise();
  Field3D particle_source(BoutReal &t);
  Field3D electron_heat_source(BoutReal &t);
  Field3D ion_heat_source(BoutReal &t);
//   int lower_source_yindex;
//   int upper_source_yindex;
  BoutReal source_length;
  BoutReal particle_amplitude;
  BoutReal electron_heat_amplitude;
  BoutReal ion_heat_amplitude;
  BoutReal on_time;
  BoutReal off_time;
  BoutReal particle_transient_amplitude;
  BoutReal electron_heat_transient_amplitude;
  BoutReal ion_heat_transient_amplitude;
private:
  bindex* position;
  int lower_source_yindex;
  int upper_source_yindex;
  bool particle_during_transient;
  bool particle_started;
  bool electron_heat_during_transient;
  bool electron_heat_started;
  bool ion_heat_during_transient;
  bool ion_heat_started;
  Field3D particle_source_field;
  Field3D electron_heat_source_field;
  Field3D ion_heat_source_field;
};

Sources::Sources() {
  position = new bindex;
  particle_during_transient = false;
  particle_started = false;
  electron_heat_during_transient = false;
  electron_heat_started = false;
  ion_heat_during_transient = false;
  ion_heat_started = false;
}

Sources::~Sources() {
  delete position;
}

void Sources::initialise() {
  // NB This assumes that the metric is CONSTANT!

  Coordinates *coord = mesh->coordinates();

  lower_source_yindex = (mesh->GlobalNy+1)/2-3 - int(source_length/2./sqrt(coord->g_22(2,2))/coord->dy(2,2)); // GlobalNy+1 in case GlobalNy is odd for some strange reason, then the source is still symmetrical and vanishes if source_length=0
  upper_source_yindex = mesh->GlobalNy/2-3 + int(source_length/2./sqrt(coord->g_22(2,2))/coord->dy(2,2));
}

Field3D Sources::particle_source(BoutReal &t) {
  if (t>on_time && t<off_time) {
    if (!particle_during_transient) {
      particle_source_field = 0.;
      start_index(position);
	do {
	  if (mesh->YGLOBAL(position->jy)>=lower_source_yindex && mesh->YGLOBAL(position->jy)<=upper_source_yindex) {
	    particle_source_field[*position] = (particle_amplitude + particle_transient_amplitude) * sin(PI*(1.*(mesh->YGLOBAL(position->jy)-lower_source_yindex))/(1.*(upper_source_yindex-lower_source_yindex)));
//           particle_source_field[*position] = (particle_amplitude + particle_transient_amplitude*PI/2*sin(PI*(t-on_time)/(off_time-on_time))) * sin(PI*(1.*(mesh->YGLOBAL(position->jy)-lower_source_yindex))/(1.*(upper_source_yindex-lower_source_yindex)));
	  }
	  else {
	    particle_source_field[*position] = 0.;
	  }
	} while (next_index3(position));
	particle_during_transient = true;
	particle_started = true;
    }
  }
  else {
    if (particle_during_transient || !particle_started) {
      particle_source_field = 0.;
      start_index(position);
	do {
	  if (mesh->YGLOBAL(position->jy)>=lower_source_yindex && mesh->YGLOBAL(position->jy)<=upper_source_yindex) {
	    particle_source_field[*position] = particle_amplitude * sin(PI*(1.*(mesh->YGLOBAL(position->jy)-lower_source_yindex))/(1.*(upper_source_yindex-lower_source_yindex)));
	  }
	  else {
	    particle_source_field[*position] = 0.;
	  }
	} while (next_index3(position));
	particle_during_transient = false;
	particle_started = true;
    }
  }
  return particle_source_field;
}

Field3D Sources::electron_heat_source(BoutReal &t) {
  if (t>on_time && t<off_time) {
    if (!electron_heat_during_transient) {
      electron_heat_source_field = 0.;
      start_index(position);
	do {
	  if (mesh->YGLOBAL(position->jy)>=lower_source_yindex && mesh->YGLOBAL(position->jy)<=upper_source_yindex) {
	    electron_heat_source_field[*position] = (electron_heat_amplitude + electron_heat_transient_amplitude) * sin(PI*(1.*(mesh->YGLOBAL(position->jy)-lower_source_yindex))/(1.*(upper_source_yindex-lower_source_yindex)));
//           electron_heat_source_field[*position] = (electron_heat_amplitude + electron_heat_transient_amplitude*PI/2*sin(PI*(t-on_time)/(off_time-on_time))) * sin(PI*(1.*(mesh->YGLOBAL(position->jy)-lower_source_yindex))/(1.*(upper_source_yindex-lower_source_yindex)));
	  }
	  else {
	    electron_heat_source_field[*position] = 0.;
	  }
	} while (next_index3(position));
	electron_heat_during_transient = true;
	electron_heat_started = true;
    }
  }
  else {
    if (electron_heat_during_transient || !electron_heat_started) {
      electron_heat_source_field = 0.;
      start_index(position);
	do {
	  if (mesh->YGLOBAL(position->jy)>=lower_source_yindex && mesh->YGLOBAL(position->jy)<=upper_source_yindex) {
	    electron_heat_source_field[*position] = electron_heat_amplitude * sin(PI*(1.*(mesh->YGLOBAL(position->jy)-lower_source_yindex))/(1.*(upper_source_yindex-lower_source_yindex)));
	  }
	  else {
	    electron_heat_source_field[*position] = 0.;
	  }
	} while (next_index3(position));
	electron_heat_during_transient = false;
	electron_heat_started = true;
    }
  }
  return electron_heat_source_field;
}

Field3D Sources::ion_heat_source(BoutReal &t) {
  if (t>on_time && t<off_time) {
    if (!ion_heat_during_transient) {
      ion_heat_source_field = 0;
      start_index(position);
	do {
	  if (mesh->YGLOBAL(position->jy)>=lower_source_yindex && mesh->YGLOBAL(position->jy)<=upper_source_yindex) {
	    ion_heat_source_field[*position] = (ion_heat_amplitude + ion_heat_transient_amplitude) * sin(PI*(1.*(mesh->YGLOBAL(position->jy)-lower_source_yindex))/(1.*(upper_source_yindex-lower_source_yindex)));
//         ion_heat_source_field[*position] = (ion_heat_amplitude + ion_heat_transient_amplitude*PI/2*sin(PI*(t-on_time)/(off_time-on_time))) * sin(PI*(1.*(mesh->YGLOBAL(position->jy)-lower_source_yindex))/(1.*(upper_source_yindex-lower_source_yindex)));
	  }
	  else {
	    ion_heat_source_field[*position] = 0.;
	  }
	} while (next_index3(position));
	ion_heat_during_transient = true;
	ion_heat_started = true;
    }
  }
  else {
    if (ion_heat_during_transient || !ion_heat_started) {
      ion_heat_source_field = 0;
      start_index(position);
	do {
	  if (mesh->YGLOBAL(position->jy)>=lower_source_yindex && mesh->YGLOBAL(position->jy)<=upper_source_yindex) {
	    ion_heat_source_field[*position] = ion_heat_amplitude * sin(PI*(1.*(mesh->YGLOBAL(position->jy)-lower_source_yindex))/(1.*(upper_source_yindex-lower_source_yindex)));
	  }
	  else {
	    ion_heat_source_field[*position] = 0.;
	  }
	} while (next_index3(position));
	ion_heat_during_transient = false;
	ion_heat_started = true;
    }
  }
  return ion_heat_source_field;
}
