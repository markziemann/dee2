module SearchPage.Main exposing (..)

import Browser.Events exposing (onClick, onKeyDown)
import Elastic as Elastic exposing (parse, serialize)
import Http as Http exposing (get)
import Json.Decode as Decode
import KeyBoardHelpers exposing (arrowDown, arrowUp, enter)
import Result
import Routes
import SearchPage.Helpers exposing (..)
import SearchPage.Types exposing (..)
import SharedTypes exposing (PaginationOffset, RemoteData(..), toWebData)
import Url.Builder as UB


search =
    Search


init : Model
init =
    { searchString = ""
    , searchMode = Strict
    , searchSuggestions = NotAsked
    , activeSuggestion = Nothing
    , suggestionsVisible = False
    , defaultPaginationOffset = PaginationOffset 20 0
    }


wrapAroundIdx : Int -> Int -> Int
wrapAroundIdx lengthSearchSuggestions idx =
    if lengthSearchSuggestions > 0 then
        modBy lengthSearchSuggestions idx

    else
        lengthSearchSuggestions


onlyData : Model -> ( Model, Cmd msg, Maybe OutMsg )
onlyData model =
    ( model, Cmd.none, Nothing )


noResults =
    Nothing


update : Msg -> Model -> ( Model, Cmd Msg, Maybe OutMsg )
update msg model =
    let
        noChange =
            ( model, Cmd.none, Nothing )

        toOutMsg =
            \paginationOffset searchResults ->
                -- OutMsg
                { searchResults = searchResults
                , searchParameters =
                    SearchParameters
                        model.searchMode
                        model.searchString
                        paginationOffset
                }
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

        Search ((SearchParameters searchMode searchString paginationOffset) as params) ->
            let
                -- Elastic.Word makes no modification to the search string
                -- https://www.elastic.co/guide/en/elasticsearch/reference/current/query-dsl-simple-query-string-query.html#_simple_query_string_syntax
                ( route, search_string ) =
                    case searchMode of
                        Strict ->
                            ( "api/simple_query_search/"
                            , Result.withDefault (Elastic.Word searchString) (parse searchString)
                                |> serialize
                            )

                        Fuzzy ->
                            ( "api/fuzzy_search/", searchString )

                url =
                    route
                        ++ (Routes.searchResultParams params
                                |> UB.toQuery
                           )

                serverQuery =
                    get
                        { url = url
                        , expect =
                            Http.expectJson
                                (GotHttpSearchResponse paginationOffset << toWebData)
                                (decodeSearchResults paginationOffset)
                        }
            in
            ( { model | searchString = searchString }
                |> hideSuggestions
            , serverQuery
            , noResults
            )

        EnterKey ->
            case ( emptyWebData model.searchSuggestions, model.activeSuggestion ) of
                ( False, Just value ) ->
                    -- Selecting a suggestion with enter takes
                    -- precedence over searching with enter
                    update (SuggestionSelected value) model

                -- recursive call should be avoided
                ( _, _ ) ->
                    update (Search <| defaultSearchParameters model) model

        -- recursive call should be avoided
        GetSearchSuggestions value ->
            -- This is some debounce on the search string to prevent spamming ElasticSearch with queries
            if model.searchString == value then
                ( model
                , get
                    { url = "api/search_as_you_type/" ++ value
                    , expect = Http.expectJson (GotSearchSuggestions << toWebData) decodeSearchSuggestions
                    }
                , noResults
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
                        ({ model | searchString = suggestion }
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
            onlyData { model | searchMode = Strict }

        FuzzySelected string ->
            onlyData { model | searchMode = Fuzzy }

        GotHttpSearchResponse paginationOffset webData ->
            ( model
            , Cmd.none
            , Just <| toOutMsg paginationOffset webData
            )

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
