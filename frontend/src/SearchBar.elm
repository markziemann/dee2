module SearchBar exposing (..)

import Array exposing (Array)
import Browser.Events exposing (onClick, onKeyDown)
import Elastic as Elastic exposing (parse, serialize)
import Http as Http exposing (get)
import Json.Decode as Decode
import KeyBoardHelpers exposing (arrowDown, arrowUp, enterKey, try)
import Keyboard.Event exposing (considerKeyboardEvent)
import Result
import SearchBarHelpers exposing (..)
import SearchBarTypes as SearchBarTypes exposing (..)
import Routes
import Nav



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
    , searchMode = Strict
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


onlyData : Model -> ( Model, Cmd msg, Maybe OutMsg )
onlyData model =
    ( model, Cmd.none, Nothing )


noResults =
    Nothing


update : Msg -> Model -> ( Model, Cmd Msg, Maybe OutMsg )
update msg model =
    let
        toOutMsg =
            \result ->
                Result.map (\{ hits, rows } -> OutMsg hits rows model.searchString) result
                    |> Result.toMaybe
    in
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
                (route, search_string) =
                    case model.searchMode of
                        Strict ->
                            ( "/simple_query_search/"
                            , Result.withDefault (Elastic.Word model.searchString) (parse model.searchString)
                                |> serialize
                            )
                        Fuzzy ->
                            ("/fuzzy_search/", model.searchString)

                serverQuery =
                    get
                        { url = route ++ search_string
                        , expect = Http.expectJson GotHttpSearchResponse decodeSearchResults
                        }
            in
            ( model |> hideSuggestions |> isWaiting, serverQuery, noResults )

        EnterKey ->
            case ( Array.isEmpty model.searchSuggestions, model.activeSuggestion ) of
                ( False, Just value ) ->
                    -- Selecting a suggestion with enter takes
                    -- precedence over searching with enter
                    update (SuggestionSelected value) model

                -- recursive call should be avoided
                ( _, _ ) ->
                    update Search model

        -- recursive call should be avoided
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
                    updateActiveSuggestion model 0
                        |> showSuggestions
                        |> onlyData

                Just value ->
                    onlyData (wrapAroundIdx model (value + 1))

        StrictSelected string ->
            onlyData {model| searchMode = Strict}

        FuzzySelected string ->
            onlyData {model| searchMode = Fuzzy}

        GotHttpSearchResponse result ->
            ( notWaiting model
            , Cmd.none
            , toOutMsg result
            )

        ClickOutOfSuggestions ->
            onlyData (hideSuggestions model)


subscriptions : Sub Msg
subscriptions =
    Sub.batch
        [ onKeyDown (considerKeyboardEvent (try [ arrowUp ArrowUp, arrowDown ArrowDown ]))
        , onClick (Decode.succeed ClickOutOfSuggestions)
        , onKeyDown (considerKeyboardEvent (enterKey EnterKey))
        ]
