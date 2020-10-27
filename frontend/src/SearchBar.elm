module SearchBar exposing (..)

import Array exposing (Array)
import Browser.Events exposing (onKeyDown)
import Elastic as Elastic exposing (parse, serialize)
import Http as Http exposing (get)
import Keyboard.Event exposing (considerKeyboardEvent)
import SearchBarHelpers
    exposing
        ( clearSearchSuggestions
        , decodeSearchSuggestions
        , delay
        , toKey
        , decodeSearchResults
        )
import SearchBarTypes as SearchBarTypes exposing (..)
import Result



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
                ( { model
                    | searchString = value
                    , activeSuggestion = Nothing
                    , searchSuggestions = Array.empty
                  }
                , Cmd.none
                )

            else
                ( { model
                    | searchString = value
                    , activeSuggestion = Nothing
                  }
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
            ( clearSearchSuggestions model, serverQuery )

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
                ( model, Cmd.none )

        GotSearchSuggestions result ->
            ( { model | searchSuggestions = Result.withDefault Array.empty result }, Cmd.none )

        SuggestionSelected value ->
            let
                searchString =
                    case Array.get value model.searchSuggestions of
                        Just string ->
                            string

                        Nothing ->
                            model.searchString
            in
            ( { model | searchSuggestions = Array.empty, searchString = searchString }, Cmd.none )

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
            ( { model | activeSuggestion = activeSuggestion }, Cmd.none )

        GotHttpSearchResponse result ->
            ({model | searchResults = Result.withDefault Array.empty result}, Cmd.none)

        ResultClicked function ->
            ({model | searchResults = function model.searchResults}, Cmd.none)



subscriptions : Model -> Sub Msg
subscriptions model =
    onKeyDown (considerKeyboardEvent toKey)
