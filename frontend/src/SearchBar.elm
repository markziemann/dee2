module SearchBar exposing (..)

import Array exposing (Array)
import Browser.Events exposing (onKeyDown)
import Elastic as Elastic exposing (parse, serialize)
import Http as Http exposing (get)
import Json.Decode as Decode
import Keyboard.Event exposing (considerKeyboardEvent)
import Result
import SearchBarHelpers exposing (..)
import SearchBarTypes as SearchBarTypes exposing (..)



-- Allows exporting of imported types


type alias Model =
    SearchBarTypes.Model


type alias Msg =
    SearchBarTypes.Msg


init : Model
init =
    { searchString = ""
    , searchSuggestions = Array.empty
    , activeSuggestion = Nothing
    , suggestionsVisible = True
    , searchResults = Array.empty
    }


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    case msg of
        SearchUpdate value ->
            if value == "" then
                -- Don't wait to clear the suggestions if there is nothing
                -- in the search box. Having search suggestions with an empty
                -- searchString causes issues for the highlightMatchingText function
                ( { model | searchString = value }
                    |> clearActiveSuggestion
                    |> clearSearchSuggestions
                , Cmd.none
                )

            else
                ( { model | searchString = value }
                    |> showSuggestions
                    |> clearActiveSuggestion
                , delay 1000 (GetSearchSuggestions value)
                )

        Search ->
            let
                -- Elastic.Word makes no modification to the search string
                search_string =
                    Result.withDefault (Elastic.Word model.searchString) (parse model.searchString)
                        |> serialize

                serverQuery =
                    get
                        { url = "/search/" ++ search_string
                        , expect = Http.expectJson GotHttpSearchResponse decodeSearchResults
                        }
            in
            ( model |> clearSearchSuggestions, serverQuery )

        GetSearchSuggestions value ->
            -- This is some debounce on the search string to prevent spamming ElasticSearch with queries
            if model.searchString == value then
                ( model
                , get
                    { url = "/search_as_you_type/" ++ value
                    , expect = Http.expectJson GotSearchSuggestions decodeSearchSuggestions
                    }
                )

            else
                ( model |> showSuggestions, Cmd.none )

        GotSearchSuggestions result ->
            ( { model | searchSuggestions = Result.withDefault Array.empty result }
                |> showSuggestions
            , Cmd.none
            )

        SuggestionSelected value ->
            let
                searchString =
                    case Array.get value model.searchSuggestions of
                        Just string ->
                            string

                        Nothing ->
                            model.searchString
            in
            ( { model | searchString = searchString }
                |> clearSearchSuggestions
            , Cmd.none
            )

        KeyPressed key ->
            let
                activeSuggestion =
                    if model.searchSuggestions /= Array.empty then
                        Just
                            (modBy (Array.length model.searchSuggestions)
                                (case model.activeSuggestion of
                                    Nothing ->
                                        case key of
                                            ArrowUp ->
                                                Array.length model.searchSuggestions

                                            ArrowDown ->
                                                0

                                    Just value ->
                                        case key of
                                            ArrowUp ->
                                                value - 1

                                            ArrowDown ->
                                                value + 1
                                )
                            )

                    else
                        Nothing
            in
            ( { model | activeSuggestion = activeSuggestion } |> showSuggestions
            , Cmd.none
            )

        GotHttpSearchResponse result ->
            ( { model | searchResults = Result.withDefault Array.empty result }, Cmd.none )

        ResultClicked function ->
            ( { model | searchResults = function model.searchResults }, Cmd.none )

        ClickOutOfSuggestions ->
            ( model |> hideSuggestions, Cmd.none )


subscriptions : Model -> Sub Msg
subscriptions model =
    Sub.batch
        [ onKeyDown (considerKeyboardEvent toKey)
        , Browser.Events.onClick (Decode.succeed ClickOutOfSuggestions)
        ]
