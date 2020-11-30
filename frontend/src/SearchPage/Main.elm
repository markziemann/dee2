module SearchPage.Main exposing (..)

import Browser.Events exposing (onClick, onKeyDown)
import Browser.Navigation as Nav
import Http as Http exposing (get)
import Json.Decode as Decode
import KeyBoardHelpers exposing (arrowDown, arrowUp, enter)
import Routes
import SearchPage.Helpers exposing (..)
import SearchPage.Types exposing (..)
import SharedTypes exposing (PaginationOffset, RemoteData(..), toWebData)
import Url.Builder


search =
    Search


init : Nav.Key -> Model
init navKey =
    let
        defaultPaginationOffset =
            PaginationOffset 20 0
    in
    { navKey = navKey
    , previousSearch = Nothing
    , level = Nothing
    , mode = Strict
    , query = ""
    , searchSuggestions = NotAsked
    , activeSuggestion = Nothing
    , suggestionsVisible = False
    , defaultPaginationOffset = defaultPaginationOffset
    }


wrapAroundIdx : Int -> Int -> Int
wrapAroundIdx lengthSearchSuggestions idx =
    if lengthSearchSuggestions > 0 then
        modBy lengthSearchSuggestions idx

    else
        lengthSearchSuggestions


onlyData : Model -> ( Model, Cmd msg )
onlyData model =
    ( model, Cmd.none )


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    let
        noChange =
            ( model, Cmd.none )
    in
    case msg of
        SearchUpdate string ->
            let
                newModel =
                    { model | query = string }
            in
            if string == "" then
                -- Don't wait to clear the suggestions if there is nothing
                -- in the search box. Having search suggestions with an empty
                -- query causes issues for the highlightMatchingText function
                onlyData
                    (newModel
                        |> clearActiveSuggestion
                        |> clearSearchSuggestions
                    )

            else
                case newModel |> toSearchParameters of
                    Just searchParameters ->
                        ( newModel |> showSuggestions
                        , delay 1000 (GetSearchSuggestions searchParameters)
                        )

                    Nothing ->
                        noChange

        Search ->
            case model |> toSearchParameters of
                Just searchParameters ->
                    ( { model | previousSearch = Just searchParameters }
                        |> hideSuggestions
                    , Nav.pushUrl model.navKey (Routes.searchResultsRoute searchParameters)
                    )

                Nothing ->
                    noChange

        EnterKey ->
            case ( emptyWebData model.searchSuggestions, model.activeSuggestion ) of
                ( False, Just value ) ->
                    -- Selecting a suggestion with enter takes
                    -- precedence over searching with enter
                    update (SuggestionSelected value) model

                -- recursive call should be avoided
                ( _, _ ) ->
                    update Search model

        -- recursive call should be avoided
        GetSearchSuggestions ((SearchParameters level mode query _) as searchParameters) ->
            -- This is some debounce on the search string to prevent spamming ElasticSearch with queries
            if query == model.query then
                ( model |> showSuggestions
                , get
                    { url =
                        "api/suggestions/"
                            ++ (Routes.searchResultParams searchParameters |> Url.Builder.toQuery)
                    , expect = Http.expectJson (GotSearchSuggestions << toWebData) decodeSearchSuggestions
                    }
                )

            else
                onlyData (showSuggestions model)

        GotSearchSuggestions webData ->
            onlyData
                ({ model | searchSuggestions = webData }
                    |> showSuggestions
                )

        SuggestionSelected index ->
            case getWebDataIndex index model.searchSuggestions of
                Just suggestion ->
                    onlyData
                        ({ model | query = suggestion }
                            |> clearActiveSuggestion
                            |> hideSuggestions
                        )

                Nothing ->
                    onlyData model

        ArrowUp ->
            case ( model.activeSuggestion, lengthWebData model.searchSuggestions ) of
                ( Nothing, Just length ) ->
                    updateActiveSuggestion model length
                        |> showSuggestions
                        |> onlyData

                ( Just value, Just length ) ->
                    updateActiveSuggestion model (wrapAroundIdx length (value - 1))
                        |> onlyData

                ( _, _ ) ->
                    noChange

        ArrowDown ->
            case ( model.activeSuggestion, lengthWebData model.searchSuggestions ) of
                ( Nothing, Just length ) ->
                    updateActiveSuggestion model 0
                        |> showSuggestions
                        |> onlyData

                ( Just value, Just length ) ->
                    updateActiveSuggestion model (wrapAroundIdx length (value + 1))
                        |> onlyData

                ( _, _ ) ->
                    noChange

        StrictSelected string ->
            onlyData { model | mode = Strict}

        FuzzySelected string ->
            onlyData { model | mode = Fuzzy}

        ClickOutOfSuggestions ->
            onlyData (hideSuggestions model)


subscriptions : Model -> Sub Msg
subscriptions model =
    [ onKeyDown <|
        Decode.oneOf
            [ enter EnterKey
            , arrowUp ArrowUp
            , arrowDown ArrowDown
            ]
    ]
        |> (\subs ->
                if model.suggestionsVisible then
                    List.append [ onClick <| Decode.succeed ClickOutOfSuggestions ] subs

                else
                    subs
           )
        |> Sub.batch
