module UI.KeyboardInputSpec exposing (main)

import Array
import Json.Encode
import Main
import MainTypes
import Process
import Runner
import SearchBarTypes
import Spec exposing (..)
import Spec.Claim as Claim
import Spec.Command
import Spec.Markup as Markup
import Spec.Markup.Event as Event
import Spec.Markup.Selector exposing (..)
import Spec.Observer as Observer
import Spec.Setup as Setup
import Spec.Step
import Spec.Time as Time
import Extra


enterKeyPressSpec : Spec.Spec MainTypes.Model MainTypes.Msg
enterKeyPressSpec =
    describe "Enter key press selects active search suggestion"
        [ Spec.scenario "Search suggestions are visibale and a selection has been activated"
            (Spec.given
                (Setup.initForApplication (wrapMock mockSearchSuggestions (Main.init ()))
                    |> Setup.withDocument Main.view
                    |> Setup.withUpdate Main.update
                    |> Setup.withSubscriptions Main.subscriptions
                )
                |> Spec.when "Search String is entered and arrow keys select suggestion and enter is pressed"
                    [

                     Markup.target << document
                    , Event.trigger "keydown" <| Json.Encode.object [ ( "key", Json.Encode.string "ArrowUp" ) ]

                    --, Spec.Command.send (Spec.Command.fake (MainTypes.GotSearchBarMsg SearchBarTypes.EnterKey))
                    ]
                |> Spec.it "selects the active suggestion"
                    (Observer.observeModel (.searchBar >> .test)
                        |> Spec.expect (Extra.equals True)
                    )
            )
        ]



--upArrowPressed =
--    Json.Encode.object [ ( "keyCode", Json.Encode.int 38 ) ]
--        |> Event.trigger "keydown"


upArrowPressed : Spec.Step.Context model -> Spec.Step.Command msg
upArrowPressed =
    Json.Encode.object [ ( "key", Json.Encode.string "ArrowUp" ) ]
        |> Event.trigger "keypress"


keyPressEvent : String -> Spec.Step.Context model -> Spec.Step.Command msg
keyPressEvent char =
    Json.Encode.object
        [ ( "key", Json.Encode.string char )
        ]
        |> Event.trigger "keypress"


enterPressed =
    Json.Encode.object [ ( "keyCode", Json.Encode.int 13 ) ]
        |> Event.trigger "keydown"


type alias ModelCmdMsg =
    ( MainTypes.Model, Cmd MainTypes.Msg )


wrapMock :
    (ModelCmdMsg -> ModelCmdMsg)
    -> (a -> b -> ModelCmdMsg)
    -> (a -> b -> ModelCmdMsg)
wrapMock func init =
    \url key ->
        func (init url key)


mockSearchSuggestions : ( MainTypes.Model, Cmd MainTypes.Msg ) -> ( MainTypes.Model, Cmd MainTypes.Msg )
mockSearchSuggestions ( model, cmd ) =
    let
        searchBar =
            model.searchBar

        newSearchBar =
            { searchBar | searchSuggestions = Array.fromList [ "1", "2", "3", "4" ] }
    in
    ( { model | searchBar = newSearchBar }, cmd )


main =
    Runner.browserProgram
        [ enterKeyPressSpec
        ]
