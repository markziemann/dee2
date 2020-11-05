module Routes exposing (..)

import SharedTypes
import Url
import Url.Builder as UB exposing (QueryParameter)
import Url.Parser as UP exposing ((</>), (<?>))
import Url.Parser.Query as Query


type Route
    = HomeRoute
    | SearchResultsRoute (Maybe String) (Maybe String) (Maybe String)


type Page
    = HomePage Route
    | SearchResultsPage Route


homeSlug =
    "/"


searchResultsSlug =
    "Search"


routeParser : UP.Parser (Route -> a) a
routeParser =
    UP.oneOf
        [ UP.map HomeRoute (UP.s homeSlug)
        , UP.map SearchResultsRoute
            (UP.s searchResultsSlug
                <?> Query.string "searchString"
                <?> Query.string "perPage"
                <?> Query.string "offset"
            )
        ]


searchResultsRoute : String -> SharedTypes.PaginationOffset -> String
searchResultsRoute searchString =
    (UB.absolute [ searchResultsSlug ]) << (searchResultParams searchString)


searchResultParams : String -> SharedTypes.PaginationOffset -> List QueryParameter
searchResultParams searchString { perPage, offset } =
    [ UB.string "searchString" searchString
    , UB.string "perPage" <| String.fromInt perPage
    , UB.string "offset" <| String.fromInt offset
    ]


homePage =
    HomePage HomeRoute


determinePage : Url.Url -> Page
determinePage url =
    case UP.parse routeParser url of
        Just page ->
            case page of
                HomeRoute ->
                    homePage

                SearchResultsRoute _ _ _ ->
                    SearchResultsPage page

        Nothing ->
            homePage
