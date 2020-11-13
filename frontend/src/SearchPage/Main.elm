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


search =
    Search


init : Nav.Key -> Model
init navKey =
    let
        defaultPaginationOffset =
            PaginationOffset 20 0
    in
    { navKey = navKey
    , searchParameters = SearchParameters Strict "" defaultPaginationOffset
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


noResults =
    Nothing


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    let
        noChange =
            ( model, Cmd.none )
    in
    case msg of
        SearchUpdate searchString ->
            let
                newSearchParameters =
                    withSearchString searchString model.searchParameters
            in
            if searchString == "" then
                -- Don't wait to clear the suggestions if there is nothing
                -- in the search box. Having search suggestions with an empty
                -- searchString causes issues for the highlightMatchingText function
                onlyData
                    ({ model | searchParameters = newSearchParameters }
                        |> clearActiveSuggestion
                        |> clearSearchSuggestions
                    )

            else
                ( { model | searchParameters = newSearchParameters }
                    |> showSuggestions
                , delay 1000 (GetSearchSuggestions newSearchParameters)
                )

        Search searchParameters ->
            ( model
                |> hideSuggestions
            , Nav.pushUrl model.navKey (Routes.searchResultsRoute searchParameters)
            )

        EnterKey ->
            case ( emptyWebData model.searchSuggestions, model.activeSuggestion ) of
                ( False, Just value ) ->
                    -- Selecting a suggestion with enter takes
                    -- precedence over searching with enter
                    update (SuggestionSelected value) model

                -- recursive call should be avoided
                ( _, _ ) ->
                    update (Search model.searchParameters) model

        -- recursive call should be avoided
        GetSearchSuggestions ((SearchParameters _ searchString _) as searchParameters) ->
            -- This is some debounce on the search string to prevent spamming ElasticSearch with queries
            if model.searchParameters == searchParameters then
                ( model
                , get
                    { url = "api/search_as_you_type/" ++ searchString
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
                        ({ model | searchParameters = withSearchString suggestion model.searchParameters }
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
            onlyData { model | searchParameters = withSearchMode Strict model.searchParameters }

        FuzzySelected string ->
            onlyData { model | searchParameters = withSearchMode Fuzzy model.searchParameters }

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
