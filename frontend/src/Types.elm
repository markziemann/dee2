module Types exposing (..)

import Browser exposing (Document)
import Browser.Navigation as Nav
import ResultsPage.Types
import SearchPage.Types exposing (SearchParameters, SearchResults)
import SharedTypes exposing (WebData)
import Url
import HomePage.Main
import Routes


type alias Model =
    { navKey : Nav.Key
    , url : Url.Url
    , homePage: HomePage.Main.Model
    , searchPage : SearchPage.Types.Model
    , resultsPage : ResultsPage.Types.Model
    , route : Routes.Route
    }


type Msg
    = GotHomePageMsg HomePage.Main.Msg
    | GotSearchRunsPageMsg SearchPage.Types.Msg
    | GotSearchProjectsPageMsg SearchPage.Types.Msg
    | GotResultsPageMsg ResultsPage.Types.Msg
    | LinkClicked Browser.UrlRequest
    | UrlChanged Url.Url
    | GotHttpSearchResponse SearchParameters (WebData SearchResults)
