module MainTypes exposing (..)

import Browser exposing (Document)
import Browser.Navigation as Nav
import SearchBar
import SearchBarTypes
import Url
import Table

type Route
    = HomeRoute
    | SearchResultsRoute (Maybe String)

type alias PageData =
    { route : Route
    , subscriptions : List (Sub Msg)
    }


type Page
    = HomePage PageData
    | SearchResultsPage PageData


type alias Model =
    { navKey : Nav.Key
    , url : Url.Url
    , searchBar : SearchBar.Model
    , page : Page
    , searchResults : SearchBarTypes.SearchResults
    , resultsTableState: Table.State
    , resultsTableQuery: String
    }

type Msg
    = GotSearchBarMsg SearchBar.Msg
    | Search -- Message of this type will be sent to the SearchBar module
    | EnterKey
    | LinkClicked Browser.UrlRequest
    | UrlChanged Url.Url
    | ResultClicked SearchBarTypes.SearchResult
    | SetResultsTableQuery String
    | SetResultsTableState Table.State


