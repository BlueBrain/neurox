int main_ring(int argc, char** argv, char** env, std::string & path){
  (void)env; // unused
  nrnmpi_init(1, &argc, &argv);
  initnrn();
  mk_mech(path.c_str());
  mk_netcvode();
  int gids[1] = {0};
  std::string filesdat = path+"files.dat";
  nrn_setup(path.c_str(), filesdat.c_str(),nrn_need_byteswap, 0);
  t = 0;
  dt = 0.025;
  double mindelay = BBS_netpar_mindelay(10.0);
  printf("mindelay = %g\n", mindelay);
  mk_spikevec_buffer(10000);
  nrn_finitialize(1, -65.0);
  BBS_netpar_solve(100.);

  return 0;
}

void modl_reg() {
	// not right place, but plays role of nrnivmodl constructed
	// mod_func.c.
}
