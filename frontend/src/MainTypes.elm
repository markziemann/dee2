module MainTypes exposing (..)

import Browser exposing (Document)
import Browser.Navigation as Nav
import SearchBar
import SearchBarTypes
import Table
import Url
import Array
import Routes
import Bytes


type alias Model =
    { navKey : Nav.Key
    , url : Url.Url
    , searchBar : SearchBar.Model
    , page : Routes.Page
    , searchHits: Maybe Int
    , searchResultRows : Maybe (Array.Array SearchBarTypes.SearchResult)
    , resultsTableState : Table.State
    , resultsTableQuery : String
    }


type Msg
    = GotSearchBarMsg SearchBar.Msg
    | LinkClicked Browser.UrlRequest
    | UrlChanged Url.Url
    | ResultClicked SearchBarTypes.SearchResult
    | SetResultsTableQuery String
    | SetResultsTableState Table.State
