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
    let
        defaultPaginationOffset =
            PaginationOffset 20 0
    in
    { searchParameters = SearchParameters Strict "" defaultPaginationOffset
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
                , searchParameters = withPagination paginationOffset model.searchParameters
                }

        updateSearchString = \searchString ->
            { model | searchParameters = withSearchString searchString model.searchParameters }
    in
    case msg of
        SearchUpdate searchString ->
            if searchString == "" then
                -- Don't wait to clear the suggestions if there is nothing
                -- in the search box. Having search suggestions with an empty
                -- searchString causes issues for the highlightMatchingText function
                onlyData
                    (updateSearchString searchString
                        |> clearActiveSuggestion
                        |> clearSearchSuggestions
                    )

            else
                ( updateSearchString searchString
                    |> showSuggestions
                , delay 1000 (GetSearchSuggestions searchString)
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
            ( updateSearchString searchString
                |> hideSuggestions
            , serverQuery
            , Just <| OutMsg Loading params
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
        GetSearchSuggestions value ->
            -- This is some debounce on the search string to prevent spamming ElasticSearch with queries
            if getSearchString model.searchParameters  == value then
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
                        (updateSearchString suggestion
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
