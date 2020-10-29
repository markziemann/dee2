port module Runner exposing
  ( program
  , browserProgram
  , skip
  , pick
  )

import Spec exposing (Message)

port elmSpecOut : Message -> Cmd msg
port elmSpecIn : (Message -> msg) -> Sub msg
port elmSpecPick : () -> Cmd msg

config : Spec.Config msg
config =
  { send = elmSpecOut
  , listen = elmSpecIn
  }

pick =
  Spec.pick elmSpecPick

skip =
  Spec.skip

program =
  Spec.program config

browserProgram =
  Spec.browserProgram config