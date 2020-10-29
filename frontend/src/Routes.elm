module Routes exposing (..)
import Url
import Url.Parser as UrlP exposing ((</>), (<?>))
import Url.Parser.Query as Query

type Route
    = HomeRoute
    | SearchResultsRoute (Maybe String)


type Page
    = HomePage Route
    | SearchResultsPage Route


homeSlug =
    "/"


searchResultsSlug =
    "SearchResults"


routeParser : UrlP.Parser (Route -> a) a
routeParser =
    UrlP.oneOf
        [ UrlP.map HomeRoute (UrlP.s homeSlug)
        , UrlP.map SearchResultsRoute (UrlP.s searchResultsSlug <?> Query.string "q")
        ]


homePage =
    HomePage HomeRoute


searchResultsPage : Maybe String -> Page
searchResultsPage maybeString =
    SearchResultsPage (SearchResultsRoute maybeString)


determinePage : Url.Url -> Page
determinePage url =
    case UrlP.parse routeParser url of
        Just page ->
            case page of
                HomeRoute ->
                    homePage

                SearchResultsRoute maybeString ->
                    searchResultsPage maybeString

        Nothing ->
            homePage
