#ifndef _RUNEVTOBJ_H_
#define _RUNEVTOBJ_H_

class runevtobj
{
  public:
    runevtobj();
    runevtobj(int run, int evt);
    ~runevtobj();
    int run_num;
    int evt_num;
    bool is_lowest;
  private:  
};
#endif
