#include <RAT/CountProc.hh>
#include <RAT/Log.hh>

namespace RAT {

CountProc::CountProc() : Processor("count") {
  dscount = 0;
  evcount = 0;
  updateInterval = 1;  // Print count message for every event
}

CountProc::~CountProc() {
  /* Do nothing */

  info << dformat("CountProc: Total # of events %d (%d triggered events)\n", dscount, evcount);
}

void CountProc::SetI(std::string param, int value) {
  if (param == "update")
    if (value > 0)
      updateInterval = value;
    else
      throw ParamInvalid(param, "update interval must be > 0");
  else
    throw ParamUnknown(param);
}

Processor::Result CountProc::DSEvent(DS::Root *ds) {
  /*
  RAT::DS::MC *mc = ds -> GetMC();
  int nb_tracks = mc -> GetMCTrackCount();
  RAT::info << "Number of tracks: " << nb_tracks << "\n";
  for (int track_id=0; track_id<nb_tracks; track_id++) {
    RAT::DS::MCTrack *track = mc -> GetMCTrack(track_id);
    int nb_steps = track -> GetMCTrackStepCount();
    RAT::info << "Number of steps: " << nb_steps << "\n";
    for (int step_id=0; step_id<nb_steps; step_id++) {
      RAT::DS::MCTrackStep *step = track -> GetMCTrackStep(step_id);
      RAT::info << "Track KE: " << step -> GetKE() << ", quenched energy: " << step -> GetScintEdepQuenched() << "\n";
    }
  }
  */
  dscount++;
  evcount += ds->GetEVCount();

  if (dscount % updateInterval == 0) info << dformat("CountProc: Event %d (%d triggered events)\n", dscount, evcount);
  return OK;

}

}  // namespace RAT
