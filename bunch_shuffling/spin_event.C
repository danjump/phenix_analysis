#include "spin_event.h" 

spin_event::spin_event() {
  filled = false;
}

spin_event::~spin_event() {
  
}

void spin_event::fill(int run, int evt, int cross, int pass_arm, int charge, int eta, int spin) {
  run_num = run;
  evt_num = evt;
  clockcross = cross;
  arm = pass_arm;
  charge_index = charge;
  eta_index = eta;
  spin_config = spin;

  filled = true;
}

void spin_event::copy(spin_event * other) {
  run_num =      other->get_run_num();
  evt_num =      other->get_evt_num();
  clockcross =   other->get_clockcross();
  arm =          other->get_arm();
  charge_index = other->get_charge_index();
  eta_index =    other->get_eta_index();
  spin_config =  other->get_spin_config();

  filled = other->is_filled();
}

bool spin_event::is_filled() {
  return filled;
}


int spin_event::get_run_num() {
  return run_num;
}

void spin_event::set_run_num(int run) {
  run_num = run;
}


int spin_event::get_evt_num() {
  return evt_num;
}

void spin_event::set_evt_num(int evt) {
  evt_num = evt;
}

int spin_event::get_clockcross() {
  return clockcross;
}

void spin_event::set_clockcross(int cross) {
  clockcross = cross;
}


int spin_event::get_arm() {
  return arm;
}

void spin_event::set_arm(int pass_arm) {
  arm = pass_arm;
}


int spin_event::get_charge_index() {
  return charge_index;
}

void spin_event::set_charge_index(int charge) {
  charge_index = charge;
}


int spin_event::get_eta_index() {
  return eta_index;
}

void spin_event::set_eta_index(int eta) {
  eta_index = eta;
}


int spin_event::get_spin_config() {
  return spin_config;
}

void spin_event::set_spin_config(int spin) {
  spin_config = spin;
}

