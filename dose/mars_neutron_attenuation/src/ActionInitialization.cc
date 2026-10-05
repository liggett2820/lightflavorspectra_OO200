#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

void ActionInitialization::Build() const {
  PrimaryGeneratorAction* primaryGenerator = new PrimaryGeneratorAction();
  SetUserAction(primaryGenerator);
  m_primaryGenerator = primaryGenerator;

  RunAction* runAction = new RunAction();
  SetUserAction(runAction);
  m_runAction = runAction;

  EventAction* eventAction = new EventAction(runAction);
  SetUserAction(eventAction);

  SetUserAction(new SteppingAction(eventAction));
}
