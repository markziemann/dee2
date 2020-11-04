module Types exposing (..)

import Browser exposing (Document)
import Browser.Navigation as Nav
import ResultsPage.Types
import SearchPage.Types
import Url
import Routes


type alias Model =
    { navKey : Nav.Key
    , url : Url.Url
    , searchPage : SearchPage.Types.Model
    , resultsPage : ResultsPage.Types.Model
    , page : Routes.Page
    }


type Msg
    = GotSearchPageMsg SearchPage.Types.Msg
    | GotResultsPageMsg ResultsPage.Types.Msg
    | LinkClicked Browser.UrlRequest
    | UrlChanged Url.Url
