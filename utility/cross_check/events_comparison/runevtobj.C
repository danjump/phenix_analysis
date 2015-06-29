#include "runevtobj.h"

runevtobj::runevtobj() {
  is_lowest = false;
}

runevtobj::runevtobj(int run, int evt) {
  is_lowest = false;
  run_num = run;
  evt_num = evt;
}

runevtobj::~runevtobj() { }


