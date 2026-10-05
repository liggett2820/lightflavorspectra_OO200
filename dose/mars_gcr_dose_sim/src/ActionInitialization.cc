#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

void ActionInitialization::Build() const {
  PrimaryGeneratorAction* primaryGenerator = new PrimaryGeneratorAction();
  SetUserAction(primaryGenerator);
  m_primaryGenerator = primaryGenerator;

  RunAction* runAction = new RunAction();
  SetUserAction(runAction);
  m_runAction = runAction;

  SetUserAction(new SteppingAction(runAction));
}
