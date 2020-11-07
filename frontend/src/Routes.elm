module Routes exposing (..)

import Maybe.Extra as ME
import SearchPage.Types exposing (SearchMode(..))
import SharedTypes
import Url
import Url.Builder as UB exposing (QueryParameter)
import Url.Parser as UP exposing ((</>), (<?>))
import Url.Parser.Query as Query


type Route
    = HomeRoute
    | SearchRoute
        { searchMode : SearchMode
        , searchString : String
        , paginationOffset : SharedTypes.PaginationOffset
        }
    | Unknown


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
            SearchRoute
                { searchMode = searchMode
                , searchString = searchString
                , paginationOffset = paginationOffset
                }

        ( _, _, _ ) ->
            Unknown


routeParser : UP.Parser (Route -> a) a
routeParser =
    UP.oneOf
        [ UP.map HomeRoute UP.top
        , UP.map parseSearchResultRoute
            (UP.s searchResultsSlug
                <?> Query.string "searchMode"
                <?> Query.string "searchString"
                <?> Query.string "perPage"
                <?> Query.string "offset"
            )
        ]


searchResultsRoute : SearchMode -> String -> SharedTypes.PaginationOffset -> String
searchResultsRoute searchMode searchString =
    UB.absolute [ searchResultsSlug ] << searchResultParams searchMode searchString


searchResultParams : SearchMode -> String -> SharedTypes.PaginationOffset -> List QueryParameter
searchResultParams searchMode searchString { perPage, offset } =
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
    , UB.string "perPage" <| String.fromInt perPage
    , UB.string "offset" <| String.fromInt offset
    ]


determinePage : Url.Url -> Route
determinePage url =
    UP.parse routeParser url
        |> ME.unwrap Unknown identity
