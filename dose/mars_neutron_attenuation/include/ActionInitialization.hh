#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class PrimaryGeneratorAction;
class RunAction;

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
