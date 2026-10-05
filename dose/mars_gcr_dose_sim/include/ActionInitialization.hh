#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class PrimaryGeneratorAction;
class RunAction;

// Serial-execution only (see RunAction.hh) -- Build() constructs the
// PrimaryGeneratorAction and RunAction and stashes pointers to them
// (Build() is const, hence "mutable") so main.cc's campaign loop can reach
// them after RunManager::Initialize() to set the current (species, energy,
// physics-variant) before each beamOn() call.
class ActionInitialization : public G4VUserActionInitialization
{
public:
  ActionInitialization() = default;
  ~ActionInitialization() override = default;

  void Build() const override;

  PrimaryGeneratorAction* GetPrimaryGeneratorAction() const { return m_primaryGenerator; }
  RunAction* GetRunAction() const { return m_runAction; }

private:
  mutable PrimaryGeneratorAction* m_primaryGenerator = nullptr;
  mutable RunAction* m_runAction = nullptr;
};

#endif
