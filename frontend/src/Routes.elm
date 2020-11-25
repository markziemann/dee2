module Routes exposing (..)

import Maybe.Extra as ME
import SearchPage.Types exposing (Level(..), Mode(..), SearchParameters(..))
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


parseLevel : Maybe String -> Maybe Level
parseLevel maybeLevel =
    case maybeLevel of
        Just "Projects" ->
            Just Projects

        Just "Runs" ->
            Just Runs

        _ ->
            Nothing


parseMode : Maybe String -> Maybe Mode
parseMode maybeSearchMode =
    case maybeSearchMode of
        Just "Strict" ->
            Just Strict

        Just "Fuzzy" ->
            Just Fuzzy

        _ ->
            Nothing


parseSearchResultRoute : Maybe String -> Maybe String -> Maybe String -> Maybe String -> Maybe String -> Route
parseSearchResultRoute maybeLevel maybeSearchMode maybeSearchString maybePerPage maybeOffset =
    case
        Maybe.map4 SearchParameters
            (parseLevel maybeLevel)
            (parseMode maybeSearchMode)
            maybeSearchString
            (Maybe.map2 SharedTypes.PaginationOffset
                (Maybe.andThen String.toInt maybePerPage)
                (Maybe.andThen String.toInt maybeOffset)
            )
    of
        Just searchParameters ->
            ResultsRoute <| searchParameters

        Nothing ->
            Unknown


routeParser : UP.Parser (Route -> a) a
routeParser =
    UP.oneOf
        [ UP.map HomeRoute <| UP.oneOf [ UP.top, UP.s "improved_search" ]
        , UP.map parseSearchResultRoute
            (UP.s searchResultsSlug
                <?> Query.string "level"
                <?> Query.string "mode"
                <?> Query.string "query"
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
searchResultParams (SearchParameters level mode query paginationOffset) =
    let
        levelString =
            case level of
                Projects ->
                    "Projects"

                Runs ->
                    "Runs"

        modeString =
            case mode of
                Strict ->
                    "Strict"

                Fuzzy ->
                    "Fuzzy"
    in
    [ UB.string "level" levelString
    , UB.string "mode" modeString
    , UB.string "query" query
    , UB.string "perPage" <| String.fromInt paginationOffset.perPage
    , UB.string "offset" <| String.fromInt paginationOffset.offset
    ]


determinePage : Url.Url -> Route
determinePage url =
    UP.parse routeParser url
        |> ME.unwrap Unknown identity
