module ResultsPage.Types exposing (..)

import Array
import Dict
import SearchPage.Types
import Set
import SharedTypes
import Table


type alias SelectedResult =
    ( Int, ( String, String ) )


type alias SelectedResults =
    Dict.Dict Int ( String, String )


type alias ResultsPendingRemoval =
    Set.Set Int


type alias ResultRows =
    Array.Array SearchPage.Types.SearchResult


type alias Model =
    { searchHits : Maybe Int
    , searchResultRows : Maybe ResultRows
    , resultsTableQuery : String
    , resultsTableState : Table.State
    , selectedResultsTableQuery : String
    , selectedResultsTableState : Table.State
    , downloading : Bool
    , selectedResults : SelectedResults
    , resultsPendingRemoval : ResultsPendingRemoval
    , paginationOffset : SharedTypes.PaginationOffset
    }


type Msg
    = ResultClicked SearchPage.Types.SearchResult
    | RemoveStagedSelections
    | SelectedResultClicked Int
    | SetResultsTableQuery String
    | SetResultsTableState Table.State
    | SetSelectedResultsTableQuery String
    | SetSelectedResultsTableState Table.State
    | DownloadRequested
    | DownloadButtonReset
    | PageRequest SharedTypes.PaginationOffset


type alias OutMsg =
    Maybe SharedTypes.PaginationOffset