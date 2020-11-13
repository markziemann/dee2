module Types exposing (..)

import Browser exposing (Document)
import Browser.Navigation as Nav
import ResultsPage.Types
import SearchPage.Types exposing (SearchParameters, SearchResults)
import SharedTypes exposing (WebData)
import Url
import Routes


type alias Model =
    { navKey : Nav.Key
    , url : Url.Url
    , searchPage : SearchPage.Types.Model
    , resultsPage : ResultsPage.Types.Model
    , route : Routes.Route
    }


type Msg
    = GotSearchPageMsg SearchPage.Types.Msg
    | GotResultsPageMsg ResultsPage.Types.Msg
    | LinkClicked Browser.UrlRequest
    | UrlChanged Url.Url
    | GotHttpSearchResponse SearchParameters (WebData SearchResults)
