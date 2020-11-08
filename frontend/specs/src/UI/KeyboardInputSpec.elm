module UI.KeyboardInputSpec exposing (main)

import Array
import Extra exposing (equals)
import Json.Encode
import Main
import Types
import Runner
import SearchPage.Types
import Spec exposing (..)
import Spec.Claim as Claim
import Spec.Markup as Markup
import Spec.Markup.Event as Event
import Spec.Markup.Selector exposing (..)
import Spec.Setup as Setup

import SharedTypes exposing (RemoteData(..))


enterKeyPressSpec : Spec.Spec Types.Model Types.Msg
enterKeyPressSpec =
    describe "Enter key press selects active search suggestion"
        [ Spec.scenario "Search suggestions are visible and a selection has been activated"
            (Spec.given
                (Setup.initForApplication (wrapMock mockSearchSuggestions (Main.init ()))
                    |> Setup.withDocument Main.view
                    |> Setup.withUpdate Main.update
                    |> Setup.withSubscriptions Main.subscriptions
                )
                |> Spec.when "Search String is entered and arrow keys select suggestion and enter is pressed"
                    [ Markup.target << by [ id "search-bar" ]
                    , Event.input "This should equal 3"
                    , Markup.target << document
                    , arrowDownPressed
                    , arrowDownPressed
                    , arrowDownPressed
                    , arrowDownPressed
                    , arrowUpPressed -- 4-1 = 3
                    , enterPressed
                    ]
                |> Spec.it "selects the active suggestion"
                    (Markup.observeElement
                        |> Markup.query << by [id "search-bar"]
                        |> Spec.expect(Claim.isSomethingWhere <| Markup.text <| equals "3")
                    )
            )
        ]


arrowUpPressed =
    keyPressed "ArrowUp"


arrowDownPressed =
    keyPressed "ArrowDown"


enterPressed =
    keyPressed "Enter"


keyPressed key =
    Json.Encode.object [ ( "key", Json.Encode.string key ) ]
        |> Event.trigger "keydown"


type alias ModelCmdMsg =
    ( Types.Model, Cmd Types.Msg )


wrapMock :
    (ModelCmdMsg -> ModelCmdMsg)
    -> (a -> b -> ModelCmdMsg)
    -> (a -> b -> ModelCmdMsg)
wrapMock func init =
    \url key ->
        func (init url key)


mockSearchSuggestions : ( Types.Model, Cmd Types.Msg ) -> ( Types.Model, Cmd Types.Msg )
mockSearchSuggestions ( model, cmd ) =
    let
        searchBar =
            model.searchPage

        newSearchBar =
            { searchBar | searchSuggestions =  Success <| Array.fromList [ "1", "2", "3", "4" ]}
    in
    ( { model | searchPage = newSearchBar }, cmd )


main =
    Runner.browserProgram
        [ enterKeyPressSpec
        ]
