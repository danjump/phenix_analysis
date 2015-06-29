#ifndef _SPIN_EVENT_H_
#define _SPIN_EVENT_H_

class spin_event
{
  public:
    spin_event();
    ~spin_event();
    
    void fill(int run, int evt, int cross, int pass_arm, int charge, int eta, int spin);
    void copy(spin_event * other);
    bool is_filled();
    int get_run_num();
    void set_run_num(int run);
    int get_evt_num();
    void set_evt_num(int evt);
    int get_clockcross();
    void set_clockcross(int cross);
    int get_arm();
    void set_arm(int pass_arm);
    int get_charge_index();
    void set_charge_index(int charge);
    int get_eta_index();
    void set_eta_index(int eta);
    int get_spin_config();
    void set_spin_config(int spin);

  private:
    bool filled;
    int run_num;
    int evt_num;
    int clockcross;
    int arm;
    int charge_index;
    int eta_index;
    int spin_config;
};


#endif //_SPIN_EVENT_H_
