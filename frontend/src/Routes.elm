module Routes exposing (..)

import Maybe.Extra as ME
import SearchPage.Types exposing (SearchMode(..), SearchParameters(..))
import SharedTypes
import Url
import Url.Builder as UB exposing (QueryParameter)
import Url.Parser as UP exposing ((</>), (<?>))
import Url.Parser.Query as Query


type Route
    = HomeRoute
    | SearchRunsRoute
    | SearchProjectsRoute
    | ResultsRoute SearchParameters
    | Unknown


searchRunsRoute =
    "Runs"

searchProjectsRoute =
    "Projects"

searchResultsSlug =
    "Search"


parseSearchMode : Maybe String -> Maybe SearchMode
parseSearchMode maybeSearchMode =
    case maybeSearchMode of
        Just "Strict" ->
            Just Strict

        Just "Fuzzy" ->
            Just Fuzzy

        _ ->
            Nothing


parseSearchResultRoute : Maybe String -> Maybe String -> Maybe String -> Maybe String -> Route
parseSearchResultRoute maybeSearchMode maybeSearchString maybePerPage maybeOffset =
    case
        ( parseSearchMode maybeSearchMode
        , maybeSearchString
        , Maybe.map2 SharedTypes.PaginationOffset
            (Maybe.andThen String.toInt maybePerPage)
            (Maybe.andThen String.toInt maybeOffset)
        )
    of
        ( Just searchMode, Just searchString, Just paginationOffset ) ->
            ResultsRoute <|
                SearchParameters searchMode searchString paginationOffset

        ( _, _, _ ) ->
            Unknown


routeParser : UP.Parser (Route -> a) a
routeParser =
    UP.oneOf
        [ UP.map HomeRoute <| UP.oneOf [ UP.top, UP.s "improved_search" ]
        , UP.map parseSearchResultRoute
            (UP.s searchResultsSlug
                <?> Query.string "searchMode"
                <?> Query.string "searchString"
                <?> Query.string "perPage"
                <?> Query.string "offset"
            )
        , UP.map SearchRunsRoute (UP.s searchRunsRoute)
        , UP.map SearchProjectsRoute (UP.s searchProjectsRoute)
        ]


searchResultsRoute : SearchParameters -> String
searchResultsRoute searchUrlParameters =
    UB.absolute [ searchResultsSlug ] <| searchResultParams searchUrlParameters


searchResultParams : SearchParameters -> List QueryParameter
searchResultParams (SearchParameters searchMode searchString paginationOffset) =
    let
        searchModeString =
            case searchMode of
                Strict ->
                    "Strict"

                Fuzzy ->
                    "Fuzzy"
    in
    [ UB.string "searchMode" searchModeString
    , UB.string "searchString" searchString
    , UB.string "perPage" <| String.fromInt paginationOffset.perPage
    , UB.string "offset" <| String.fromInt paginationOffset.offset
    ]


determinePage : Url.Url -> Route
determinePage url =
    UP.parse routeParser url
        |> ME.unwrap Unknown identity
