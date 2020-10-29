module SearchBar exposing (..)

import Array exposing (Array)
import Browser.Events exposing (onKeyDown)
import Elastic as Elastic exposing (parse, serialize)
import Http as Http exposing (get)
import Json.Decode as Decode
import KeyBoardHelpers exposing (arrowDown, arrowUp, try)
import Keyboard.Event exposing (considerKeyboardEvent)
import Maybe.Extra as MExtra
import Result
import SearchBarHelpers exposing (..)
import SearchBarTypes as SearchBarTypes exposing (..)



-- Allows exporting of imported types


type alias Model =
    SearchBarTypes.Model


type alias Msg =
    SearchBarTypes.Msg



-- Expose this Msg constructor to other modules


searchMsg =
    SearchBarTypes.Search


suggestionSelectedMsg =
    SearchBarTypes.SuggestionSelected


init : Model
init =
    { searchString = ""
    , searchSuggestions = Array.empty
    , activeSuggestion = Nothing
    , suggestionsVisible = True
    , waitingForResponse = False
    }


wrapAroundIdx : Model -> Int -> Model
wrapAroundIdx model idx =
    let
        arrayLength =
            Array.length model.searchSuggestions
    in
    if arrayLength > 0 then
        modBy arrayLength idx
            |> updateActiveSuggestion model

    else
        model


onlyData : Model -> ( Model, Cmd msg, Maybe SearchResults )
onlyData model =
    ( model, Cmd.none, Nothing )


noResults =
    Nothing


update : Msg -> Model -> ( Model, Cmd Msg, Maybe SearchResults )
update msg model =
    case msg of
        SearchUpdate value ->
            if value == "" then
                -- Don't wait to clear the suggestions if there is nothing
                -- in the search box. Having search suggestions with an empty
                -- searchString causes issues for the highlightMatchingText function
                onlyData
                    ({ model | searchString = value }
                        |> clearActiveSuggestion
                        |> clearSearchSuggestions
                    )

            else
                ( { model | searchString = value }
                    |> showSuggestions
                , delay 1000 (GetSearchSuggestions value)
                , noResults
                )

        Search ->
            let
                -- Elastic.Word makes no modification to the search string
                -- https://www.elastic.co/guide/en/elasticsearch/reference/current/query-dsl-simple-query-string-query.html#_simple_query_string_syntax
                search_string =
                    Result.withDefault (Elastic.Word model.searchString) (parse model.searchString)
                        |> serialize

                serverQuery =
                    get
                        { url = "/search/" ++ search_string
                        , expect = Http.expectJson GotHttpSearchResponse decodeSearchResults
                        }
            in
            ( model |> hideSuggestions |> isWaiting, serverQuery, noResults )

        GetSearchSuggestions value ->
            -- This is some debounce on the search string to prevent spamming ElasticSearch with queries
            if model.searchString == value then
                ( model
                , get
                    { url = "/search_as_you_type/" ++ value
                    , expect = Http.expectJson GotSearchSuggestions decodeSearchSuggestions
                    }
                , noResults
                )

            else
                onlyData (showSuggestions model)

        GotSearchSuggestions result ->
            onlyData
                ({ model | searchSuggestions = Result.withDefault Array.empty result }
                    |> showSuggestions
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
            onlyData
                ({ model | searchString = searchString }
                    |> clearActiveSuggestion
                    |> hideSuggestions
                )

        ArrowUp ->
            --onlyData model
            case model.activeSuggestion of
                Nothing ->
                    updateActiveSuggestion model (Array.length model.searchSuggestions)
                        |> showSuggestions
                        |> onlyData

                Just value ->
                    onlyData (wrapAroundIdx model (value - 1))

        ArrowDown ->
            case model.activeSuggestion of
                Nothing ->
                    (updateActiveSuggestion model 0)
                        |> showSuggestions
                        |> onlyData

                Just value ->
                    onlyData (wrapAroundIdx model (value + 1))

        GotHttpSearchResponse result ->
            let
                searchResults =
                    Result.withDefault Array.empty result
            in
            ( notWaiting model
            , Cmd.none
            , Just searchResults
            )

        ClickOutOfSuggestions ->
            onlyData (hideSuggestions model)


subscriptions : Sub Msg
subscriptions =
    Sub.batch
        [ onKeyDown (considerKeyboardEvent (try [ arrowUp ArrowUp, arrowDown ArrowDown ])) ]



--[ onKeyDown (considerKeyboardEvent (arrowDown ArrowDown))
--, Browser.Events.onClick (Decode.succeed ClickOutOfSuggestions)
--]
--[]
