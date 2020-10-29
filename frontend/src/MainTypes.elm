module MainTypes exposing (..)

import Browser exposing (Document)
import Browser.Navigation as Nav
import SearchBar
import SearchBarTypes
import Table
import Url


type Route
    = HomeRoute
    | SearchResultsRoute (Maybe String)


type Page
    = HomePage Route
    | SearchResultsPage Route


type alias Model =
    { navKey : Nav.Key
    , url : Url.Url
    , searchBar : SearchBar.Model
    , page : Page
    , searchResults : SearchBarTypes.SearchResults
    , resultsTableState : Table.State
    , resultsTableQuery : String
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
